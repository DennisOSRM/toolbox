mod command_line;
use std::{fs::File, io::BufWriter};

use bincode::serialize_into;
use env_logger::Env;
use log::info;

use crate::command_line::{Arguments, InputFormat};
use toolbox_rs::{ddsg, dimacs, edge::InputEdge, metis};

fn main() {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    println!(r#"     ___     _       _                    "#);
    println!(r#"    | _ \   | |     (_)     ___      _ _  "#);
    println!(r#"    |  _/   | |     | |    / -_)    | '_| "#);
    println!(r#"   _|_|_   _|_|_   _|_|_   \___|   _|_|_  "#);
    println!(r#" _| """ |_|"""""|_|"""""|_|"""""|_|"""""| "#);
    println!(r#" "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-' "#);

    // parse and print command line parameters
    let args = <Arguments as clap::Parser>::parse();
    info!("{args}");

    let edges: Vec<InputEdge<i32>> = match args.input_format {
        InputFormat::DDSG => ddsg::read_graph(&args.graph, ddsg::WeightType::Original),
        InputFormat::DIMACS => dimacs::read_graph(&args.graph, dimacs::WeightType::Original),
        InputFormat::METIS => metis::read_graph(&args.graph, metis::WeightType::Original),
    };

    let coordinates = match args.input_format {
        InputFormat::DDSG => ddsg::read_coordinates(&args.coordinates),
        InputFormat::DIMACS => dimacs::read_coordinates(&args.coordinates),
        InputFormat::METIS => metis::read_coordinates(&args.coordinates),
    };

    info!("writing edges into intermediate format");
    let mut f = BufWriter::new(File::create(args.graph + ".toolbox").unwrap());
    serialize_into(&mut f, &edges).unwrap();

    info!("writing coordinates into intermediate format");
    let mut f = BufWriter::new(File::create(args.coordinates + ".toolbox").unwrap());
    serialize_into(&mut f, &coordinates).unwrap();

    info!("done.");
}
