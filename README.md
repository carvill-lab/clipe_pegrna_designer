<img src="src/www/clipe_logo.png" alt="CliPE Logo" height="75">

## Overview

CliPE pegRNA designer is a tool designed to facilitate the construction of pegRNA sequences for targeted prime editing screens. This tool simplifies the process of designing pegRNAs by automating the sequence generation and validation steps.

## Features

- Automated prioritization of editing regions based on Clinvar and GnomAD data
- Automated pegRNA sequence generation
- Sequence validation and error checking

## Usage

CliPE pegRNA designer is available for use at [design.clipe-mave.org](https://design.clipe-mave.org). Additional information on usage can be found on the CliPE homepage [here](https://home.clipe-mave.org) or in the manuscript [here]().

If you have used our design tools for your prime editing experiments, please cite our manuscript as follows: ____

### Running Locally
If you would like to run CliPE pegRNA designer locally, you can pull our docker image from Github and run with the following instructions. This should take no longer than 5 minutes to setup

1. Download and install [Docker Desktop](https://www.docker.com/get-started/), then open the application
2. Using the command line, pull our docker image: 
```
docker pull --platform linux/x86_64 ghcr.io/nbodkin/clipe_pegrna_designer:main
```
3. Using the command line, run the image:
```
docker run -p 8080:8080 --name clipe ghcr.io/nbodkin/clipe_pegrna_designer:main
```
4. Open up the webapp in your browser using http://localhost:8080

## Contributing

We welcome contributions to our web application. If you would like to contribute, please fork the repository and submit a pull request.

## License

This project is licensed under the Creative Commons Attribution 4.0 International Public License. Learn more [here](https://creativecommons.org/licenses/by/4.0/).

## Contact

For any questions or issues, please open an issue on the [GitHub repository](https://github.com/nbodkin/clipe_pegrna_designer/issues).
