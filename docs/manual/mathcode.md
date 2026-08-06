# Math Code

*Creates pieces of code to be interpreted by Giraffe to evaluate a time-series, usually employed for time-varying load/displacement prescription.*

### Syntax:

Begin code data End

-   Begin: keyword to mark the beginning of a piece of code

-   End: keyword to mark the end of a piece of code

-   code data: code to be used to evaluate the desired time-series (as a function of the variable “t” - time)

### Example:

//piece of code 1

Begin sin(2\*pi\*t) End

// piece of code 2

Begin exp(t^2) End

//piece of code 3

Begin if(t&gt;=0.5,100\*t,0) End

### Additional information:

MathCode provides Giraffe the possibility of employing the result of a piece of code (as a function of a variable “t” (time) to evaluate a desired time-series value. The code may be a simple expression, but may have a more complex structure, including multiple lines, conditional statements, loops, etc.

The syntax to be employed for coding is found in the documentation of the library “exprtk”, used inside Giraffe to interpret the pieces of code:

<http://www.partow.net/programming/exprtk/index.html>
