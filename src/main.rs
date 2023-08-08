use anyhow::Result;
use simplelog::SimpleLogger;

fn main() -> Result<()> {

    let _ = SimpleLogger::init(log::LevelFilter::Info, simplelog::Config::default());

    let gth = dform::PseudoGTH::from_cp2k_db(
        "/home/ahmad/software/cp2k/cp2k-2023.1/data/GTH_POTENTIALS",
        "Na",
        "pade",
        1
        )?;

    println!("{:#?}", gth);

    gth.to_psp8("gth.psp8", "gth diagonal form")?;

    Ok(())
}

