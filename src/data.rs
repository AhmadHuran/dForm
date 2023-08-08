use log::warn;

pub static ATOMIC_SYMBOLS: [&str; 118] = [
    "h",                                                                                                                                                                                    "he",
    "li", "be",                                                                                                                                               "b",  "c",  "n",  "o",  "f",  "ne",
    "na", "mg",                                                                                                                                               "al", "si", "p",  "s",  "cl", "ar",
    "k",  "ca",                                                                                   "sc", "ti", "v",  "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr",
    "rb", "sr",                                                                                   "y",  "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i",  "xe",
    "cs", "ba", "la","ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb","dy", "ho", "er", "tm", "yb", "lu", "hf", "ta", "w",  "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po", "at", "rn",
    "fr", "ra", "ac","th", "pa", "u",  "np", "pu", "am", "cm", "bk","cf", "es", "fm", "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds", "rg", "cn", "nh", "fl", "mc", "lv", "ts", "og",
];


pub fn atomic_number(symbol: &str) -> Option<u32> {

    let ss = symbol
        .to_lowercase()
        .trim()
        .to_owned();

    for (atom_num, atom_sym) in ATOMIC_SYMBOLS
        .iter()
        .enumerate() {
            if atom_sym == &ss {
                return Some(atom_num as u32 + 1)
            }
    }

    None

}

pub fn pspxc(functional: &str) -> i32 {
    match functional.trim().to_lowercase().as_str() {
        "pade" | "lda" => 1,
        "pbe" => 11,
        "blyp" => 18,
        "bp" => 19,
        "hcth120" => 17,
        "hcth407" => 27,
        "olyp" => 25,
        "pbesol" => 0,
        _ => {
            warn!("Unknown xc functional! \
                  PBE xc functional is assumed.");
            11
        },
    }
}
