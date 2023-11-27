# Euriklis  QUADPROG package.

This package implements the [A. Idnani and D. Goldfarb dual algorithm](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.521.6352&rep=rep1&type=pdf). It is well known that this algorithm was initially implemented from [Alberto Santini](https://www.npmjs.com/package/quadprog). We provide an exact translation of the Fortran code in a javascript - like pattern. Thus we avoid the goto_* functions which increase the time efficiency and the readability of the algorithm implementation. 
# Usage

The quadratic programming is widely used in the optimization procedures, the portfolio management with restrictions, the least squares method with restrictions and many other cases. We use this algorithm this algorithm for the so called supported vector machines (or SVM - neuronal networks architectures) - an architecture of neuronal networks. 

# Installation of the package

To install the package run in terminal:
```sh
npm install @euriklis/quadprog@latest --save
```
To use the package you have to declare the package with the conventional require approach:

```js
const solveQP = require('@euriklis/quadprog');
let Dmat = [], dvec = [], Amat = [], bvec = [], res;

Dmat[0] = [];
Dmat[1] = [];
Dmat[2] = [];
Dmat[0][0] = 1;
Dmat[1][0] = 0;
Dmat[2][0] = 0;
Dmat[0][1] = 0;
Dmat[1][1] = 1;
Dmat[2][1] = 0;
Dmat[0][2] = 0;
Dmat[1][2] = 0;
Dmat[2][2] = 1;

dvec[0] = 0;
dvec[1] = 5;
dvec[2] = 0;

Amat[0] = [];
Amat[1] = [];
Amat[2] = [];
Amat[0][0] = -4;
Amat[1][0] = -3;
Amat[2][0] = 0;
Amat[0][1] = 2;
Amat[1][1] = 1;
Amat[2][1] = 0;
Amat[0][2] = 0;
Amat[1][2] = -2;
Amat[2][2] = 1;

bvec[0] = -8;
bvec[1] = 2;
bvec[2] = 0;
console.log(solveQP(Dmat, dvec, Amat, bvec));
```

# Dependencies

The package does not use any dependencies and is implemented without using of classes, methods and additional packages/functions. This pattern implementation was preferred because of the simplicity and time efficiency, which is crucial for this class of algorithms. 

# Technical details

The implemented from [Alberto Santini package](https://github.com/albertosantini/node-quadprog) is accurately translated. However, in the original Fortran subroutine also the lagrangian multipliers are provided, so in our implementation we obtain these multipliers and returns then in the output of the function. 

# Tests

Tests of the time efficiency and the accuracy of the algorithm can be run with running of the following commands in the terminal:
```sh
git clone https://github.com/VelislavKarastoychev/tests-for-quadprog
cd tests-for-quadprog && npm t
```

# Bugs and tips

Everyone who wants to inform me about useful things and practices can sends me an email to exel_mmm@abv.bg or euriklis@hotmail.com.

# License

MIT License. This package will be provided for free to any user that use it for personal and non commercial usage. The author of the package is not liable for any errors in third party software, libraries, packages and source code used at these libraries. The author also may not responsible for some possible bugs that may exists in the library.
