BEGIN checklist of SPECS file parameters and their default values
* Printing
  Major print level      1   * one-line major iteration log
  Minor print level       0   * one-line minor iteration log
  Print file         9   *
  Summary file         6   * typically the screen
  Print frequency       100   * minor iterations log on PRINT file
  Summary frequency       100   * minor iterations log on SUMMARY file
  Solution         Yes   * on the PRINT file
* Suppress options listing       * default: options are listed
  System information       No   * prints more system information
 
* QP subproblems
  QPSolver              CG * default
  Crash option         3   * first basis is essentially triangular
  Elastic mode         No   * until it seems necessary
  Elastic weight       1.0e+4   * used only during elastic mode
  Iterations limit       10000   * or 20m if that is more
  Partial price       1   * 10 for large LPs
 
* SQP method
  Major iterations limit     10000000   * or m if that is more
  Minor iterations limit     250   *
* Major step limit       2.0   *
  Superbasics limit       n+1
  Hessian dimension       750   * or Superbasics limit if that is less
  Derivative level       3   *
  Derivative linesearch       *
 
End of SPECS file checklist 
