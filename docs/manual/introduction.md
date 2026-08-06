# Introduction

Giraffe is the acronym of “Generic Interface Readily Accessible for Finite Elements”. It is a platform coded using C++ language, with the objective of generating a base-interface to be used by researchers, to implement their own finite element formulations. Giraffe does not have the mission of being a completely generic platform, which would be too ambitious. Structural problems, however, which may include translational and rotational degrees of freedom, such as possible multiphysics applications, can be sketched in such a way that permits creating a platform to embrace new elements, new contact formulations, new constraint equations, among other features. With that aim, “Giraffe Project” was started on 2014 by Prof. Alfredo Gay Neto, at University of São Paulo, Brazil.

Giraffe has started as a generalization of a previous-developed finite element code, named “FemCable”, which had the objective of simulating offshore structures: risers for oil exploitation. It had implementations of geometric nonlinear beam elements and classical node to surface contact formulation. Since a natural expansion required including new contact models, new structural elements and other resources, Giraffe was designed to have all the models included in “FemCable”. Furthermore, it was thought to embrace easy inclusion of new resources, using object orientation programming. Giraffe is under continuous development by Prof. Alfredo Gay Neto and co-workers.

On 2018 Giraffe was completely re-structured to encompass a new broad of resources. The new code structure provides new possibilities for modeling, with a higher versatility. This includes the possibility of defining a sequence of solutions, possibly mixing static and dynamic methods, according to convenience. Furthermore, “BoolTable" keyword has replaced “Steps” keyword, no longer available. “BoolTable” provides an easier way to define in which solution step each resource will be included or not included in simulation. With that, one may straightforwardly switch on/off boundary conditions, loads, joints, contacts, etc. This leads to the possibility of creating scenarios where load sequence is an issue. Furthermore, it provides numerical strategies to achieve solution of challenging nonlinear problems. Also, post-processing possibilities have changed, with a more organized set of post files, for post-processing using Paraview<sup>TM</sup> environment. With respect to the user input file, there is a slight change between Giraffe 1.0 and Giraffe 2.0 syntax, since some commands have changed, there are new ones and others were discontinued. Users are invited to have a look at new tutorials, in order to get used to new resources.

This user’s manual has some brief explanations on how to use Giraffe to simulate structural models using the available elements, contact algorithms and special constraints. The focus is to explain the syntax of the Giraffe input file. Examples may be found in Giraffe tutorials documentation.

Suggestions are always welcome and can be emailed to Prof. Alfredo Gay Neto: <alfredo.gay@usp.br>.

Enjoy!

Alfredo Gay Neto

São Paulo, Brazil, 2018.
