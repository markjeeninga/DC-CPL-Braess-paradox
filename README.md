# DC-CPL-Braess-paradox

This code considers all test cases in MATPOWER 7.0, available at https://matpower.org.
A DC power grid is derived from each test case by considering the decoupled reactive power flow.
This means that we take the branch inductances in each case as the line resistances in the DC grid.
In this analogy, PQ buses correspond to constant-power loads and PV buses correspond to voltages sources.

The goal of this code is to answer the following question:
For which test cases does there exist a connected component in the induced subgraph of PQ buses with at least 
two nodes such that the Laplacian submatrix that describes the inductive connections between these PQ buses and the 
neighboring PV buses has full column rank?
If so, then any increase in the edges within this subgraph could violate the power flow feasibility for
some choice of power demands at the loads, and hence Braess' paradox for power flow feasibilty can occur.
