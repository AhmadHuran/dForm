use std::f64::consts::PI;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::str::Lines;

use anyhow::{Result, anyhow};

use log::warn;
use ndarray::Array1;
use ndarray::Array2;
use ndarray::Axis;
use ndarray::s;
use special::Error;

use crate::pspxc;
use crate::ELogGamma;
use crate::atomic_number;
use crate::inverse_3x3;

#[derive(Debug)]
#[allow(unused)]
pub struct PseudoGTH {
    atom_nr: u32,
    n_elec: u32,
    elec_conf: Array1<u32>,
    r_loc: f64,
    alpha: f64,
    n_gau: usize,
    c_i: Array1<f64>,
    l_max: usize,
    r_l: Array1<f64>,
    n_proj: Array1<usize>,
    h: Vec<Array2<f64>>,
    pspxc: String,
}

impl PseudoGTH {
    pub fn from_single(filename: &str) -> Result<Self> {
        let f = match File::open(filename){
            Ok(ok) => ok,
            Err(err) => {
                let mut message = err.to_string();
                message.push_str(format!("\nFile: '{}'", filename).as_str());
                return Err(anyhow!(message))
            }
        };
        let mut reader = BufReader::new(f);

        let mut lines = String::new();
        reader.read_to_string(&mut lines)?;
        Self::from_lines(lines.lines())
    }

    pub fn from_cp2k_db(
        filename: &str,
        element: &str,
        functional: &str,
        core_charge: usize,
        ) -> Result<Self> {

        let f = match File::open(filename){
            Ok(ok) => ok,
            Err(err) => {
                let mut message = err.to_string();
                message.push_str(format!("\nFile: '{}'", filename).as_str());
                return Err(anyhow!(message))
            }
        };

        let mut reader = BufReader::new(f);

        let mut lines = String::new();
        reader.read_to_string(&mut lines)?;

        let lines = lines.lines();

        let cond = |v: String| -> bool {
            let p = !v.contains(element.to_lowercase().as_str())
            || !v.contains(functional.to_lowercase().as_str())
            || !v.contains(format!("q{}", core_charge).as_str());
            p
        };
        let mut potential_string = format!(
            "{} GTH-{}-q{}\n",
            element,
            functional,
            core_charge
            );

        let lines = lines
            .skip_while(
                |line| cond(line.to_lowercase())
                )
            .skip(1)
            .take_while(
                |line| !line.to_lowercase().contains("gth")
                )
            .map(|v| {
                let mut v = v.to_owned();
                v.push('\n');
                v
            });

        potential_string.extend(lines);
        if potential_string
            .lines()
            .skip(1)
            .collect::<Vec<_>>()
            .is_empty() {
            return Err(
                anyhow!(
                    "Could not find the specified potential!\n  \
                    - element: '{}',\n  \
                    - xc functional: '{}',\n  \
                    - zion: {}\n\
                    potential database file: '{}'",
                    element,
                    functional,
                    core_charge,
                    filename,
                    )
                )
        };

        Self::from_lines(potential_string.lines())
    }

    fn from_lines(lines: Lines) -> Result<Self> {

        let mut lines = lines.filter(
            |x|
            (! x.is_empty()) ||
            (! x.starts_with('#'))
            );
        let mut header = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|v| v.to_owned())
            .collect::<Vec<_>>();
        let symbol = header.remove(0);
        let zz = match atomic_number(symbol.as_str()) {
            Some(number) => number,
            None => return Err(anyhow!("Invalid chemical symbol '{}' \
                                       while reading pseudo potential!",
                                       symbol)
                               ),
        };

        let pspxc = header
            .into_iter()
            .find(|v|
                    v.to_owned()
                    .to_lowercase()
                    .starts_with("gth-")
                    )
            .unwrap_or_else(
                || {
                    warn!("Could not parse xc functional! \
                          PBE xc functional is assumed.");
                    "x-pbe".to_owned()
                }
                )
            .split('-')
            .nth(1)
            .unwrap_or_else(
                || {
                    warn!("Could not parse xc functional! \
                          PBE xc functional is assumed.");
                    "pbe"
                }
                )
            .to_owned();

        let elec: Vec<_> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|x| match x.parse::<u32>(){
                Ok(number) => number,
                Err(_error) => panic!("Cant convert elec_conf while reading pseudo potential!"),
            }).collect();

        let elec: Array1<u32> = Array1::from_vec(elec);
        let n_elec = elec.sum();

        let mut local: Vec<_> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(String::from)
            .collect();
        let r_loc: f64 = local.remove(0).parse().unwrap();
        let n_gau: usize = local.remove(0).parse().unwrap();
        let mut c_i: Array1<f64> = Array1::from_vec(local
            .iter()
            .map(|x| x.parse::<f64>().unwrap())
            .collect());
        if n_gau == 0 {
            c_i = Array1::<f64>::zeros(1);
        }

        let mut l_max = lines
            .next()
            .unwrap()
            .split_whitespace()
            .take(1)
            .map(|v| v.parse::<usize>()
                .unwrap_or_else(
                    |_| {
                        warn!("Could not parse lmax! \
                              lmax = 0 is assumed.");
                        1
                    })
                )
            .next()
            .unwrap();

        let mut h = Vec::<Array2<f64>>::with_capacity(l_max + 1);
        let mut r_l = Array1::<f64>::zeros(l_max + 1);
        let mut n_proj  = Array1::<usize>::zeros(l_max + 1);
        for ll in 0..l_max {
            let mut h_ll = Vec::<f64>::new();
            let mut at_hand: Vec<_> = lines
                .next()
                .unwrap()
                .split_whitespace()
                .map(String::from)
                .collect();
            r_l[ll] = at_hand.remove(0).parse::<f64>().unwrap();
            let n_proj_ll = at_hand.remove(0).parse::<usize>().unwrap();
            for h1j in at_hand.iter(){
                h_ll.push(h1j.parse::<f64>().unwrap());
            }
            for _proj_i in 1..n_proj_ll {
                let at_hand: Vec<_> = lines
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .map(String::from)
                    .collect();

                for h1j in at_hand.iter(){
                    h_ll.push(h1j.parse::<f64>().unwrap());
                }
            }

            let mut mat_h_ll = Array2::<f64>::zeros((n_proj_ll, n_proj_ll));

            for ii in 0..n_proj_ll {
                for jj in ii..n_proj_ll {
                    mat_h_ll[[ii, jj]] = h_ll.remove(0);
                    mat_h_ll[[jj, ii]] = mat_h_ll[[ii, jj]];
                }
            }

            n_proj[ll] = n_proj_ll;
            h.push(mat_h_ll);
        }
        let alpha = 1.0 / (2.0_f64.sqrt() * r_loc);

        if l_max == 0 {
            h.push(Array2::<f64>::zeros((1, 1)));
        }
        else {
            l_max -= 1;
        }


        Ok(
            PseudoGTH {
                atom_nr: zz.to_owned(),
                n_elec,
                elec_conf: elec,
                r_loc,
                alpha,
                n_gau,
                c_i,
                l_max,
                r_l,
                n_proj,
                h,
                pspxc,

            }
        )
    }

    pub fn hl(&self, l: usize) -> Option<&Array2<f64>> {
        if l <= self.l_max {
        return Some(&self.h[l])
        };
        None
    }

    pub fn rl(&self, l: usize) -> Option<f64> {
        if l <= self.l_max {
        return Some(self.r_l[l])
        };
        None
    }

    pub fn dump_projectors(&self, l: usize, r0: f64, r1: f64, nr: usize) -> Option<Array2<f64>> {

        let r = uniform_grid(r0, r1, nr);
        let nproj = self.n_proj[l];
        if nproj == 0 {
            return None
        }
        let rl = self.r_l[l];
        let mut dump = Array2::<f64>::zeros((nr, nproj + 1));
        dump.slice_mut(s![.., 0]).assign(&r);

        for ii in 1..=nproj {
            let fac = norm_fac(l, ii, rl);
            dump.slice_mut(s![.., ii])
                .assign(
                    &r
                    .mapv(
                        |ri| ri.powi((l + 2*ii - 1) as i32)//additional r factor because abinit
                        * (-0.5 * (ri /rl).powi(2)).exp()
                        * fac
                          )
                      );
        }

        Some(dump)
    }

    pub fn diagonal_form(&self, l: usize) -> Result<(Array2<f64>, Vec<f64>)> {
        let h = &self.h[l];
        let size = self.n_proj[l];
        if size == 0 {
            return Ok((Array2::zeros((1,1)), vec![0.0]))
        }
        let (atomic, overlap) = self.atomic_wavefuncs(l)?;
        let mut ci = Vec::<f64>::with_capacity(size);
        let mut new_wave_set = Array2::<f64>::zeros((size, size));
        new_wave_set.slice_mut(s![.., 0]).assign(&atomic.slice(s![.., 0]));
        let shs = overlap.dot(h).dot(&overlap);

        let mut new_v_old = new_wave_set.t().dot(&shs).dot(&atomic);
        let mut new_v_new = new_wave_set.t().dot(&shs).dot(&new_wave_set);

        ci.push(1.0/new_v_new[[0, 0]]);
        for ii in 1..size {
            let columns = new_wave_set.clone();
            new_wave_set.slice_mut(s![.., ii]).assign(
                &(
                    &atomic.slice(s![.., ii])
                    - columns
                    .columns()
                    .into_iter()
                    .enumerate()
                    .take(ii)
                    .map(|(jj, jcol)| &jcol * (new_v_old[[jj, ii]]*ci[jj]))
                    .reduce(|acc,v| acc + v)
                    .unwrap()
                    )
                );

        new_v_old = new_wave_set.t().dot(&shs).dot(&atomic);
        new_v_new = new_wave_set.t().dot(&shs).dot(&new_wave_set);
        ci.push(1.0/new_v_new[[ii, ii]]);

        };

        let new_proj_set = h.dot(&overlap).dot(&new_wave_set);
        Ok((new_proj_set, ci))
    }

    pub fn atomic_wavefuncs(
        &self,
        l: usize,
        ) -> Result<(Array2<f64>, Array2<f64>)> {
        let size = self.n_proj[l];
        let overlap = overlap_mat(l, self.r_l[l], size);
        let h_new = overlap
            .dot(&self.h[l])
            .dot(&overlap);
        let mut x = Array2::<f64>::zeros((size, size));
        let h_inv = inverse_3x3(&h_new).unwrap();
        for ii in 0..size {
            let xx = h_inv.dot(&overlap.slice(s![.., ii]));
            let fac = 1.0;//xx.dot(&overlap).dot(&xx).sqrt();
            x.slice_mut(s![.., ii])
                .assign(&(xx/fac));
        };

        Ok((x, overlap))
    }

    pub fn dump_diagonal_form(&self, l: usize, r0: f64, r1:f64, nr: usize) -> Result<Option<Array2<f64>>> {

        let Some(mut dump) = self.dump_projectors(l, r0, r1, nr) else {
            return Ok(None);
        };

        let (coeff, _) = self.diagonal_form(l)?;
        let projectors = dump.slice(s![.., 1..]).to_owned();
        dump.slice_mut(s![.., 1..]).assign(&(projectors.dot(&coeff)));
        Ok(Some(dump))

    }

    pub fn dump_local(&self, dr: f64, eps: f64) -> Array2<f64> {
        let r0 = 0.0;
        let r1 = self.find_bracket(eps);
        let rmax = bisection(
            (r0, r1),
            |r| {
                ((self.alpha * r).compl_error()
                - self.sr_local(r) * r / self.n_elec as f64
                ).abs() - eps
            },
            1.0e-14);
        let size = (rmax/dr + 1.0).ceil() as usize;
        let mut out = Array2::<f64>::zeros((size, 2));
        let r = uniform_grid(r0, rmax, size);
        out.slice_mut(s![.., 0]).assign(&r);
        out.slice_mut(s![.., 1]).assign(&r.mapv(|ri| -(self.n_elec as f64) / ri * (self.alpha*ri).error() + self.sr_local(ri)));
        out[[0,1]] = -(self.n_elec as f64) * 2.0 * self.alpha / PI.sqrt() + self.c_i[0];
        out

    }

    pub fn find_bracket(&self, eps: f64)  -> f64 {
        let mut r = 1.0;
        let mut val = 1.0;

        while val > 0.0 {
            r *= 2.0;
            val = (
                (self.alpha * r).compl_error()
                - self.sr_local(r) * r / self.n_elec as f64
                ).abs() - eps;

        }

        r
    }

    pub fn sr_local(&self, r: f64) -> f64 {
        if self.n_gau == 0 {
            return 0.0
        }
        let alpha = self.alpha;
        self.c_i
            .iter()
            .enumerate()
            .map(
                |(ii, ci)|
                ci
                * (2.0_f64.sqrt() * alpha * r).powi((2*ii) as i32)
                * (-(alpha*r).powi(2)).exp()
            )
            .sum::<f64>()
    }

    pub fn to_psp8(
        &self,
        file_path: &str,
        title: &str,
        ) -> Result<()> {

        let local = self.dump_local(0.01, 1.0e-15);
        let mmax = local.len_of(Axis(0));
        let rmax = local[[mmax-1, 0]];
        let mut file =  File::options()
             .write(true)
             .create_new(true)
             .open(file_path)?;

        writeln!(file, "{}", title)?;

        writeln!(
            file,
            "{:.4}\t{:.4}\t000000\t zatom,zion,pspd",
            self.atom_nr as f64,
            self.n_elec as f64
            )?;

        writeln!(
            file,
            "8\t{}\t{}\t6\t{}\t0\tpspcod,pspxc,lmax,lloc,mmax,r2well",
            pspxc(self.pspxc.as_str()),
            self.l_max,
            mmax,
            )?;

        writeln!(
            file,
            "0.0\t0.0\t0.0\trchrg,fchrg,qchrg",
            )?;

        let mut nproj = vec!["0".to_owned();5];
        self.n_proj
            .iter()
            .enumerate()
            .for_each(|(ii, n)| {nproj[ii] = n.to_string()});

        writeln!(
            file,
            "{}\tnproj",
            nproj.join("\t"),
            )?;

        writeln!(
            file,
            "0\t\t\t\textension_switch",
            )?;

        for l in 0..=self.l_max {

            let Some(dump) = self.dump_diagonal_form(l, 0.0, rmax, mmax)? else {
                continue
            };
            let (_, ekbl) = self.diagonal_form(l)?;

            writeln!(
                file,
                "{}\t\t{}",
                l,
                ekbl.into_iter().map(|v| format!("{:.15e}", v)).collect::<Vec<_>>().join("\t"),
                )?;

            for (ii, row) in dump.rows().into_iter().enumerate() {
             writeln!(
                file,
                "{}\t{}",
                ii+1,
                row.into_iter().map(|v| format!("{:.15e}", v)).collect::<Vec<_>>().join("\t"),
                )?;
            }
        }

        writeln!(file, "6")?;
            for (ii, row) in local.rows().into_iter().enumerate() {
             writeln!(
                file,
                "{}\t{}",
                ii+1,
                row.into_iter().map(|v| format!("{:.15e}", v)).collect::<Vec<_>>().join("\t"),
                )?;
            }

        Ok(())
    }
}

// i starts from 1
pub fn norm_fac(l: usize, i: usize, rl: f64) -> f64 {
    let ilh = (2 * i + l) as f64 - 0.5;
    (2.0/ilh.elgamma()).sqrt() / rl.powf(ilh)
}

pub fn overlap(l: usize, i: usize, j: usize, rl: f64) -> f64 {
    let ijlh = (i + j + l) as f64 - 0.5;
    0.5 * rl.powf(2.0*ijlh) * ijlh.elgamma()
        * norm_fac(l, i, rl)
        * norm_fac(l, j, rl)
}

pub fn overlap_mat(l: usize, rl: f64, size: usize) -> Array2<f64> {
    let mut out = Array2::<f64>::eye(size);
    for i in 0..size {
        for j in i+1..size {
            let ii = i + 1 ;
            let jj = j + 1 ;
            out[[i, j]] = overlap(l, ii, jj, rl);
            out[[j, i]] = out[[i, j]];
        }
    }
    out
}

pub fn uniform_grid(start: f64, end:f64, size:usize) -> Array1<f64> {
    let dx = (end-start)/(size as f64 - 1.0);
    Array1::from_iter((0..size).map(|i| start + i as f64 * dx))
}

pub fn bisection<F>(bracket: (f64, f64), function: F, eps: f64) -> f64
where F: Fn(f64) -> f64 {
    assert!(eps > 0.0);
    let mut fr = function(bracket.1);
    let mut fl = function(bracket.0);

    if fl * fr >= 0.0 {
        panic!("Bad bracket in bisection!")
    }

    let mut left = bracket.0;
    let mut right = bracket.1;
    let mut fm = 1.0_f64;
    let mut mid = 0.0_f64;
    while fm.abs() > eps {
        mid = (left + right)/2.0;
        fm = function(mid);

        if fm * fr < 0.0 {
            fl = fm;
            left = mid;
            continue;
        };

        if fm * fl < 0.0 {
            fr = fm;
            right = mid;
            continue;
        };

    }

    mid
}
