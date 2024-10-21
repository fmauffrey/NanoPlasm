import sys
import os
import yaml
import shutil

def check_files(folder):
    """ Verify if fastq files are present in the specified folder """
    try:
        files_list = os.listdir(folder)
        if not files_list:
            print("Specified folder is empty")
            sys.exit(2)
    except FileNotFoundError:
        print("Specified folder does not exist")
        sys.exit(2)
    
    return True


def create_config(nano_dir, ill_dir, run_dir):
    """ Create a config file for the run """

    # Load the configuration file template
    config_template = yaml.safe_load(open(f"{os.path.dirname(sys.argv[0])}/data/config_template.yaml", "r"))

    # Define reads folders path and samples names and add the info to the config file
    samples = [fastq.replace(".fastq", "") for fastq in os.listdir(nano_dir)]
    samples_nano_path = os.path.abspath(nano_dir)
    samples_ill_path = os.path.abspath(ill_dir) if ill_dir else ""
    config_file = {"path_nano": samples_nano_path, "path_ill": samples_ill_path, **config_template, "samples": samples}

    # Save new configuration file
    try:
        os.mkdir(run_dir)
    except FileExistsError:
        if input(f"Directory {run_dir} already exists. Overwrite ? [y/N] -> ") == "y":
            shutil.rmtree(run_dir)
            os.mkdir(run_dir)
        else:
            sys.exit(0)

    with open(f"{run_dir}/config.yaml", "w") as output:
        yaml.safe_dump(config_file, output, sort_keys=False)


def update_config(nano_dir, ill_dir, run_dir):
    """ Update the configuration file """

    # Load the configuration file template
    config_file = yaml.safe_load(open(f"{run_dir}/config.yaml", "r"))

    # Define reads folders path and samples names and add the info to the config file
    samples = [fastq.replace(".fastq", "") for fastq in os.listdir(nano_dir)]
    samples_nano_path = os.path.abspath(nano_dir)
    samples_ill_path = os.path.abspath(ill_dir) if ill_dir else ""
    config_file["path_nano"] = samples_nano_path
    config_file["path_ill"] = samples_ill_path
    config_file["samples"] = samples

    # Save new configuration file
    with open(f"{run_dir}/config.yaml", "w") as output:
        yaml.safe_dump(config_file, output, sort_keys=False)