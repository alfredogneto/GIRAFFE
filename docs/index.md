# Giraffe

<p align="center">
  <img src="assets/giraffe-logo.png" alt="Giraffe" width="260">
</p>

**Generic Interface Readily Accessible for Finite Elements**

Giraffe is a C++ platform for implementing and running finite element formulations, aimed at structural problems with translational and rotational degrees of freedom (and possible multiphysics applications). It has been developed since 2014 by Prof. Alfredo Gay Neto and co-workers at the University of São Paulo, Brazil.

Giraffe was originally a Windows-only code. The **macOS version** (this repository's *Mac-version* branch) was contributed by **Dr. Ajay B Harish** (Harish Research Lab, The University of Manchester) — see [Installing Giraffe](manual/general-information.md#installing-giraffe) for platform-specific build instructions. For questions about the Mac version, please [open an issue on GitHub](https://github.com/Harish-Research-Lab/GIRAFFE/issues).

This site is the browsable version of the **Giraffe User's Manual v2.0.124** and the **Giraffe Tutorials 2.0** guide, converted from the original Word/PDF documents so they are easier to search and link to. The original PDFs remain available as authoritative fallbacks: [User's Manual](assets/Giraffe-Users-Manual-v2.0.124.pdf) · [Tutorials](assets/Giraffe-Tutorials-2.0-v21.pdf).

## Getting started

- [Introduction](manual/introduction.md) — what Giraffe is and its background
- [General information](manual/general-information.md) — installing and running Giraffe, input file structure, output files
- [Tutorials](tutorials/index.md) — 16 worked examples, from a first cantilever beam to contact and mooring-line dynamics, each with a ready-to-run input file
- [Nodes](manual/nodes.md) and [Elements](manual/elements.md) — the basic building blocks of a model

## Reference

Use the navigation on the left to browse the full input-file keyword reference — organized the same way as the manual: model definition (nodes, elements, materials, sections, contact...), loads and constraints, geometry primitives, monitoring/output, and solver options.

## Building Giraffe itself

For build instructions on macOS, Linux, or Windows, see the [README](https://github.com/Harish-Research-Lab/GIRAFFE/blob/Mac-version/README.md) in the repository.
