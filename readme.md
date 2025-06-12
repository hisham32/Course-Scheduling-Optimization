This repository contains a mixed integer linear programming (MILP) optimization model to schedule courses in the Department of Civil and Environmental Engineering at the University of Utah. The primary objective is to assign courses to available time slots while respecting prerequisites, corequisites, and other constraints such as slot preferences and lab-specific requirements.

# Installation requirements
1. Register for gurobi WLS license
2. Install python 3.12
3. Run pip install -r requirements.txt in the command prompt to install all the libraries

# Running the program
1. Run generate_preference.bat to generate the preference.xlsx.
2. Run the scheduling_optimization.bat to optimize the schedule; the outputs will be in [semester]_schedule.xlsx
