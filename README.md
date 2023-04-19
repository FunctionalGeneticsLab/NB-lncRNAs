[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![DOI](https://zenodo.org/badge/517749799.svg)](https://zenodo.org/badge/latestdoi/517749799)

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/FunctionalGeneticsLab/NB-lncRNAs">
    <!-- <img src="images/pipelineSimplels.png" alt="Logo" width=400>-->
  </a>

  <h3 align="center">NB-lncRNAs</h3>

  <p align="center">
    Normal Breast lncRNAs Pipeline
    <br />
    <a href="https://github.com/FunctionalGeneticsLab/NB-lncRNAs"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <!--<a href="https://github.com/FunctionalGeneticsLab/NB-lncRNAs">View Demo</a>
    ·--->
    <a href="https://github.com/FunctionalGeneticsLab/NB-lncRNAs/issues">Report Bug</a>
    ·
    <a href="https://github.com/FunctionalGeneticsLab/NB-lncRNAs/issues">Request Feature</a>
  </p>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#overview">Overview</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#pipeline-prerequisites">Pipeline Prerequisites</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## Overview

This repository contains the scripts used for the discovery of thousands of long noncoding RNAs (lncRNAs) expressed in human mammary epithelial cells and for the analysis of their expression in individual cells of the normal breast epithelium.

Each script file foccus on one part of the analysis, from De novo assembly to the remapping of the lncRNAs found on TCGA data. Scripts 1 to 11 should be executed in this order:

1. Script File 1: "***de novo***" transcriptome pipeline.
2. Script File 2: NB-lncRNA annotation
3. Script File 3: Expression of NB-lncRNAs in the bulk RNAseq samples.
4. Script File 4: Transcript counts for single-cell clustering.
5. Script File 5: Single-cell clustering for Fluidigm and 10x Genomics.
6. Script File 6: Cluster-specificity index (CSI).
7. Script File 7: SCENT SR value calculation.
8. Script File 8: Monocle cell entropy.
9. Script File 9: Slingshot trajectories.
10. Script File 10: TCGA expression.
11. Script File 11: MGFR.
12. Script File 12: Accessory scripts.

<!-- GETTING STARTED -->
## Pipeline Prerequisites

To get a local copy up and running, make sure you have each script file prerequisites instaled and up to date.

<!--  To get a local copy up and running, make sure you have the following prerequisites instaled and up to date (according to each script file mentioned before):

1. Script File 1 -->

<!-- USAGE EXAMPLES -->
## Usage

Each script generates a logoutput file and can be executed with 
```sh
bash script-name.sh
```

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request


<!-- ACKNOWLEDGEMENTS
## Acknowledgements

* []()
* []()
* []() -->


<!-- CONTACT -->
## Contact

Isabela Almeida - mb.isabela42@gmail.com

Project Link: [https://github.com/FunctionalGeneticsLab/NB-lncRNAs](https://github.com/FunctionalGeneticsLab/NB-lncRNAs)

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[contributors-url]: https://github.com/FunctionalGeneticsLab/NB-lncRNAs/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[forks-url]: https://github.com/FunctionalGeneticsLab/NB-lncRNAs/network/members
[stars-shield]: https://img.shields.io/github/stars/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[stars-url]: https://github.com/FunctionalGeneticsLab/NB-lncRNAs/stargazers
[issues-shield]: https://img.shields.io/github/issues/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[issues-url]: https://github.com/FunctionalGeneticsLab/NB-lncRNAs/issues
[license-shield]: https://img.shields.io/github/license/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[license-url]: https://github.com/FunctionalGeneticsLab/NB-lncRNAs/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
