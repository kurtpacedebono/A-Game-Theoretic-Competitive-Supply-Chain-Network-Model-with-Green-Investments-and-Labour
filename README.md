# A Game Theoretic Competitive Supply Chain Network Model with Green Investments and Labour
By Kurt Pace Debono, Maria Kontorinaki, Monique Sciortino 

The constraints and initial solutions are stored for each scenario in the 'inputscenario_' file, where _ is replaced with the corresponding number of the scenario. 

In the 'scenario_' files, the VI equations are derived by differentiating the utility functions with respect to each variable x_p and v_l. 

To solve a scenario, the following command is entered into the command window:

     [xstar xtmat]=extragradient('inputscenario_',zeta);
where:

-xstar stores the final solution  

-xtmat stores the result at each iteration of the extragradient algorithm, which will allow analysis on how the solutions evolved     with each iteration

-inputscenario_ refers to the scenario being solved

-zeta is the step size 
