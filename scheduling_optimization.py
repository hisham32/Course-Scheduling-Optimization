# -*- coding: utf-8 -*-
import os
import re
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import itertools
import gurobipy as gp
from gurobipy import Model, GRB, quicksum

# Create an environment with your WLS license
def extract_license_parameters(lic_file_path):
    params = {}
    
    with open(lic_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Skip empty lines or comments (assuming comments start with '#')
            if not line or line.startswith('#'):
                continue
            
            # Check for key=value pairs
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()

                # Convert LICENSEID to integer if needed
                if key.upper() == 'LICENSEID':
                    try:
                        value = int(value)
                    except ValueError:
                        pass

                # Store key-value pair in the dictionary, using uppercase keys for consistency
                params[key.upper()] = value

    return params

license_path = "gurobi.lic"

if os.path.exists(license_path):
    params = extract_license_parameters(license_path)
    env = gp.Env(params=params)
    print("Gurobi environment started with license parameters.")
else:
    env = gp.Env()  # Will look for license automatically or use default behavior
    print("Gurobi environment started without explicit license parameters.")

# Prompt user for semester input
user_input = input("Enter semester (Fall or Spring): ").strip()

# Normalize input to lowercase for comparison and validate
if user_input.lower() == 'fall':
    plan_semester = 'Fall'
elif user_input.lower() == 'spring':
    plan_semester = 'Spring'
else:
    raise ValueError("Invalid input. Please enter 'Fall' or 'Spring'.")

# Define the path to the master Excel file
master_file = 'input_data.xlsx'

# Load each DataFrame from its corresponding sheet
catalog_df      = pd.read_excel(master_file, sheet_name='catalog')
flowchart_df    = pd.read_excel(master_file, sheet_name='flowchart')
timeslots50_df  = pd.read_excel(master_file, sheet_name='timeslots50')
timeslots80_df  = pd.read_excel(master_file, sheet_name='timeslots80')
df_prerequisites= pd.read_excel(master_file, sheet_name='prerequisites')
df_corequisites = pd.read_excel(master_file, sheet_name='corequisites')
df_capstone     = pd.read_excel(master_file, sheet_name='capstone')
df_design       = pd.read_excel(master_file, sheet_name='design')
df_tech         = pd.read_excel(master_file, sheet_name='tech')
df_discussion   = pd.read_excel(master_file, sheet_name='discussion', index_col=0)
df_lab          = pd.read_excel(master_file, sheet_name='lab', index_col=0)

# Filter courses offered in the plan semester
catalog_df['course'] = catalog_df['Subject code'] + ' ' + catalog_df['Catalog Number'].astype(str)
courses_offered = []
for index, row in catalog_df.iterrows():
    course = row['course']
    if course == 'CVEEN 5920':
        continue
    if row['Semesters Typically Offered']==plan_semester or row['Semesters Typically Offered']=='All':
        courses_offered.append(course)


# courses_offered = courses_offered + list(lab_courses.keys())
df_preference = pd.read_excel('preference.xlsx')

# Convert to dict for further processing
slots_50_min = timeslots50_df.to_dict(orient='list')
slots_80_min = timeslots80_df.to_dict(orient='list')
flowchart = flowchart_df.to_dict(orient='list')

# Convert slots into DataFrames for 50-minute and 80-minute
slots_50 = {'Day': [], 'Start': [], 'End': []}
days_50 = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday']
for day in days_50:
    start_key = f'{day}_start'
    end_key = f'{day}_end'
    if start_key in slots_50_min:
        start_times = [t for t in slots_50_min[start_key] if pd.notna(t)]
        end_times = [t for t in slots_50_min[end_key] if pd.notna(t)]
        slots_50['Day'].extend([day]*len(start_times))
        slots_50['Start'].extend(start_times)
        slots_50['End'].extend(end_times)
df_slots_50 = pd.DataFrame(slots_50)

slots_80 = {'Day': [], 'Start': [], 'End': []}
days_80 = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday']
for day in days_80:
    start_key = f'{day}_start'
    end_key = f'{day}_end'
    if start_key in slots_80_min:
        start_times = [t for t in slots_80_min[start_key] if pd.notna(t)]
        end_times = [t for t in slots_80_min[end_key] if pd.notna(t)]
        slots_80['Day'].extend([day]*len(start_times))
        slots_80['Start'].extend(start_times)
        slots_80['End'].extend(end_times)
df_slots_80 = pd.DataFrame(slots_80)

# Combine all time slots into one dataframe
df_slots_50['SlotLength'] = 50
df_slots_80['SlotLength'] = 80
df_slots = pd.concat([df_slots_50, df_slots_80], ignore_index=True)
df_slots['SlotID'] = df_slots.index

def time_to_minutes(t):
    h, m, s = map(int, t.split(':'))
    return h * 60 + m

df_slots['StartMin'] = df_slots['Start'].apply(time_to_minutes)
df_slots['EndMin'] = df_slots['End'].apply(time_to_minutes)

# Initialize model
model = gp.Model("Course Scheduling",env=env)
# model.params.MIPFocus = 1  # Prioritize finding good feasible solutions
# model.params.Heuristics = 0.5  # Prioritize finding good feasible solutions
# model.setParam('PoolSearchMode', 2)     # Aggressive search for more solutions
# model.setParam('PoolSolutions', 3)      # Request up to 3 solutions in the solution pool

# Sets
slots = df_slots['SlotID'].tolist()
days = df_slots['Day'].unique().tolist()

day_numbers = {'Monday': 1, 'Tuesday': 2, 'Wednesday': 3, 'Thursday': 4, 'Friday': 5}
day_short = {'Monday':'Mo', 'Tuesday':'Tu', 'Wednesday':'We', 'Thursday':'Th', 'Friday':'Fr'}

slot_lengths = df_slots.set_index('SlotID')['SlotLength'].to_dict()
slot_days = df_slots.set_index('SlotID')['Day'].to_dict()
slot_start_times = df_slots.set_index('SlotID')['Start'].to_dict()
slot_start_mins = df_slots.set_index('SlotID')['StartMin'].to_dict()
slot_end_mins = df_slots.set_index('SlotID')['EndMin'].to_dict()

capstone_courses = df_capstone['course'].tolist()
design_courses = df_design['course'].tolist()
tech_electives = df_tech['course'].tolist()

lab_courses_dict = df_lab.to_dict(orient='dict')

# Filter out only valid courses
courses = list(set([c for c in courses_offered if (c not in lab_courses_dict) and (1000 <= int(c.split()[1]) < 5999)]))
# 3000-level courses
courses_3000 = [c for c in courses_offered if (c not in lab_courses_dict) and (3000 <= int(c.split()[1]) < 4000)]

#########################################
# Build Multi-Section Lab Courses
#########################################
lab_sections_map = {}  # { 'CVEEN 3015': ['CVEEN 3015_sec1', 'CVEEN 3015_sec2'], ... }
lab_courses = []
for main_lab, num_secs in lab_courses_dict.items():
    section_names = []
    for i in range(1, num_secs['sections'] + 1):
        sec_name = f"{main_lab}_sec{i}"
        section_names.append(sec_name)
        courses.append(sec_name)  # Add it to the global courses list
        lab_courses.append(sec_name)
        courses_3000.append(sec_name)
    lab_sections_map[main_lab] = section_names

#########################################
# Build Multi-Section Discussion Courses
#########################################
def get_discussion_course(course):
    # Return the discussion course name if it exists, else None
    for main_c, (disc_c, ctype, sections) in required_discussions.items():
        if main_c == course:
            return [f"{disc_c}_sec{i}" for i in range(1, sections + 1)]
    return None

def create_discussion_line(course_list):
    """
    Return multiple "variations" of the original course_list,
    replacing ONLY the first course that has discussion with each of its
    discussion sections. Ignore subsequent courses that also have discussions.
    """
    discussion_lines = []
    replaced = False

    for i, c in enumerate(course_list):
        # Check if c has discussion sections
        dcourse = get_discussion_course(c)

        # Only replace the first course that has discussion
        if not replaced and dcourse is not None and all(sec in courses for sec in dcourse):
            # Create a new variation for each discussion section
            for disc_sec in dcourse:
                new_line = course_list.copy()
                new_line[i] = disc_sec
                discussion_lines.append(new_line)
            replaced = True  # Mark that we've done a replacement

        # If we've already replaced a course, do nothing further for subsequent courses

    # If we never replaced anything (no course had a valid discussion),
    # just return the original list as one entry
    if not replaced:
        discussion_lines.append(course_list)

    return discussion_lines

# Identify main courses that need a discussion
required_discussions = df_discussion.to_dict(orient='list')

discussion_to_main = {}
for main_c, (disc_c, ctype, sections) in required_discussions.items():
    for i in range(1, sections + 1):
        discussion_to_main[f"{disc_c}_sec{i}"] = main_c

# Add discussion courses only if their main course is offered
discussion_courses = {}
for main_course, (disc_course, ctype, sections) in required_discussions.items():
    if main_course in courses_offered:
        for i in range(1, sections + 1):
            disc_course_sec = f"{disc_course}_sec{i}"
            # The main course is offered, so we offer the discussion course
            discussion_courses[disc_course_sec] = {
                'type': ctype,
                '50_min_needed': 3 if ctype == 'mixed' else (1 if ctype == '50_only' else 0),
                '80_min_needed': 2 if ctype == 'mixed' else (1 if ctype == '80_only' else 0),
            }

discussion_courses_map = {}  # store detail about each discussion main + all sections
for main_c, (disc_c, disc_type, disc_num_secs) in required_discussions.items():
    if main_c in courses_offered:
        # Create multiple sections for disc_c
        sec_names = []
        for i in range(1, disc_num_secs + 1):
            sec_name = f"{disc_c}_sec{i}"
            sec_names.append(sec_name)
            courses.append(sec_name)
        # Store the type + the # of 50/80 needed for each section
        discussion_courses_map[disc_c] = {
            'type': disc_type,
            'num_secs': disc_num_secs,
            'sections': sec_names
        }

x = model.addVars(courses, slots, vtype=GRB.BINARY, name='x')
y50 = model.addVars(courses, vtype=GRB.BINARY, name='y50')
y80 = model.addVars(courses, vtype=GRB.BINARY, name='y80')
z = model.addVars(courses, vtype=GRB.BINARY, name='z')

# Identify distinct start times for slot lengths
distinct_start_times_50 = sorted(set(df_slots.loc[df_slots['SlotLength'] == 50, 'StartMin']))
distinct_start_times_80 = sorted(set(df_slots.loc[df_slots['SlotLength'] == 80, 'StartMin']))

slots_for_start_50 = {t: [] for t in distinct_start_times_50}
slots_for_start_80 = {t: [] for t in distinct_start_times_80}

for s in slots:
    length = slot_lengths[s]
    start_min = slot_start_mins[s]
    if length == 50:
        slots_for_start_50[start_min].append(s)
    else:
        slots_for_start_80[start_min].append(s)

# Overlaps
overlaps = {}
for s1 in slots:
    for s2 in slots:
        if s1 >= s2:
            continue
        day1 = slot_days[s1]
        day2 = slot_days[s2]
        if day1 != day2:
            continue
        start1 = slot_start_mins[s1]
        end1 = slot_end_mins[s1]
        start2 = slot_start_mins[s2]
        end2 = slot_end_mins[s2]
        if max(start1, start2) < min(end1, end2):
            overlaps.setdefault(s1, set()).add(s2)
            overlaps.setdefault(s2, set()).add(s1)

# consider for discussion
# Prerequisite synchronization
u = {}
for index, row in df_prerequisites.iterrows():
    prereq_courses = [row[str(i)] for i in range(1, len(row)) if pd.notna(row[str(i)])]
    prereq_courses = [c for c in prereq_courses if c in courses]
    if len(prereq_courses) <= 1:
        continue
    for s in slots:
        u_vars = [x[c, s] for c in prereq_courses]
        u_name = f"Prereq_{'_'.join(prereq_courses)}_{s}"
        u_var = model.addVar(vtype=GRB.BINARY, name=u_name)
        model.addGenConstrMin(u_var, u_vars)
        u[(tuple(prereq_courses), s)] = u_var

    # Discussion line
    discussion_lines = create_discussion_line(prereq_courses)
    for discussion_line in discussion_lines:
        # If discussion_line differs from prereq_courses (or even if not), apply constraints:
        if len(discussion_line) > 1:  # Still meaningful if at least 2 courses
            for s in slots:
                u_vars = [x[c, s] for c in discussion_line]
                u_name = f"Prereq_{'_'.join(discussion_line)}_{s}"
                u_var = model.addVar(vtype=GRB.BINARY, name=u_name)
                model.addGenConstrMin(u_var, u_vars)
                u[(tuple(discussion_line), s)] = u_var

# Early slot penalty variables
early_slots = [s for s in slots if slot_start_times[s] == '07:30:00']

# Friday-late slots penaly variables
friday_late_slots = []
for s in slots:
    if slot_days[s] == "Friday":
        start_min = slot_start_mins[s]
        if start_min >= 14*60:  # 14:00 = 14*60 minutes
            friday_late_slots.append(s)
penalty_slots = early_slots + friday_late_slots
# Create penalty variables for 07:30:00 or "friday after 14:00" usage
w_penalty = model.addVars(courses, penalty_slots, vtype=GRB.BINARY, name="wFriLate")


w_overlap_3000 = {}
w_overlap_design = {}
w_overlap_tech = {}

for course in courses:
    if course in discussion_courses.keys():
        continue
    # Day variables for each course
    d = model.addVars(days, vtype=GRB.BINARY, name=f"d_{course}")

    # Slot length choice
    # Course preference
    if (course in lab_courses) or (course not in df_preference['Course'].tolist()):
        pref_50 = 1
        pref_80 = 1
    else:
        pref_50 = df_preference.loc[df_preference['Course'] == course, '50_min'].values[0]
        pref_80 = df_preference.loc[df_preference['Course'] == course, '80_min'].values[0]


    if course in lab_courses:
        model.addConstr(z[course] == 0, name=f"Preference_{course}_lab")
        model.addConstr(y50[course] + y80[course] == 1, name=f"SlotLengthChoice_{course}")
    elif course in capstone_courses:
        model.addConstr(z[course] == 0, name=f"Preference_{course}_lab")
        model.addConstr(y50[course] == 0, name=f"{course}_No50")
        model.addConstr(y80[course] == 1, name=f"{course}_Forced80")
    else:
        model.addConstr(y50[course] + y80[course] == 1, name=f"SlotLengthChoice_{course}")
        if pref_50 == 1 and pref_80 == 0:
            model.addConstr(z[course] <= y50[course], name=f"Preference_{course}_50")
        elif pref_80 == 1 and pref_50 == 0:
            model.addConstr(z[course] <= y80[course], name=f"Preference_{course}_80")
        else:
            model.addConstr(z[course] <= y50[course] + y80[course], name=f"Preference_{course}_both")

    # If lab course
    if course in lab_courses:
        n_slots_50 = 3
        n_slots_80 = 2
        sequence_vars = []
        sequence_info = []
        # Identify consecutive sequences for lab courses
        for day in days:
            # 50-min sequences
            slots_50_day = [ss for ss in slots if slot_days[ss] == day and slot_lengths[ss] == 50]
            slots_50_day_sorted = sorted(slots_50_day, key=lambda ss: slot_start_mins[ss])
            for idx in range(len(slots_50_day_sorted) - n_slots_50 + 1):
                seq_slots = slots_50_day_sorted[idx:idx + n_slots_50]
                times = [slot_start_mins[ss] for ss in seq_slots]
                if all(times[i+1] - times[i] == 55 for i in range(len(times)-1)):
                    seq_var = model.addVar(vtype=GRB.BINARY, name=f"seq_{course}_{day}_50_{idx}")
                    sequence_vars.append(seq_var)
                    sequence_info.append({'day': day, 'length': 50, 'slots': seq_slots, 'var': seq_var})

            # 80-min sequences
            slots_80_day = [ss for ss in slots if slot_days[ss] == day and slot_lengths[ss] == 80]
            slots_80_day_sorted = sorted(slots_80_day, key=lambda ss: slot_start_mins[ss])
            for idx in range(len(slots_80_day_sorted) - n_slots_80 + 1):
                seq_slots = slots_80_day_sorted[idx:idx + n_slots_80]
                times = [slot_start_mins[ss] for ss in seq_slots]
                if all(times[i+1] - times[i] == 95 for i in range(len(times)-1)):
                    seq_var = model.addVar(vtype=GRB.BINARY, name=f"seq_{course}_{day}_80_{idx}")
                    sequence_vars.append(seq_var)
                    sequence_info.append({'day': day, 'length': 80, 'slots': seq_slots, 'var': seq_var})

        # One sequence chosen
        model.addConstr(gp.quicksum(sequence_vars) == 1, name=f"LabCourse_Sequence_{course}")

        # Link sequences to x
        for seq in sequence_info:
            seq_var = seq['var']
            seq_slots = seq['slots']
            for ss in seq_slots:
                model.addConstr(x[course, ss] == seq_var, name=f"SeqAssign_{course}_{ss}")

        # Total slots match chosen sequence
        model.addConstr(
            gp.quicksum(x[course, ss] for ss in slots) == gp.quicksum(len(seq['slots']) * seq['var'] for seq in sequence_info),
            name=f"LabCourse_TotalSlots_{course}"
        )

        # Link day variables
        for day in days:
            sequences_on_day = [seq['var'] for seq in sequence_info if seq['day'] == day]
            model.addConstr(
                d[day] == gp.quicksum(sequences_on_day),
                name=f"LabCourse_Day_{course}_{day}"
            )
        model.addConstr(
            gp.quicksum(d[day] for day in days) == 1,
            name=f"LabCourse_SingleDay_{course}"
        )

    else:
        if course == "CVEEN 2000":
            # Force CVEEN 2000 to use only one 50-min slot
            model.addConstr(y50[course] == 1, name=f"{course}_Forced50")
            model.addConstr(y80[course] == 0, name=f"{course}_No80")
            # Also, for a non-lab, the default code sets #slots=3*y50 => 3 for 50-min,
            # so override it with a custom constraint that it meets exactly once:
            model.addConstr(
                gp.quicksum(x[course, s] for s in slots if slot_lengths[s] == 50) == 1,
                name=f"{course}_One50Slot"
            )
        else:
            # Non-lab courses
            model.addConstr(
                gp.quicksum(x[course, ss] for ss in slots if slot_lengths[ss] == 50) == 3 * y50[course],
                name=f"SlotCount50_{course}"
            )
            model.addConstr(
                gp.quicksum(x[course, ss] for ss in slots if slot_lengths[ss] == 80) == 2 * y80[course],
                name=f"SlotCount80_{course}"
            )
            # Once per day
            for day in days:
                model.addConstr(
                    gp.quicksum(x[course, ss] for ss in slots if slot_days[ss] == day) <= 1,
                    name=f"OncePerDay_{course}_{day}"
                )

            d[course] = model.addVars(days, vtype=GRB.BINARY, name=f"d_{course}")
            for day in days:
                slots_on_day = [ss for ss in slots if slot_days[ss] == day]
                model.addConstr(
                    d[course][day] == gp.quicksum(x[course, ss] for ss in slots_on_day),
                    name=f"DayAssignment_{course}_{day}"
                )

            model.addConstr(
                gp.quicksum(d[course][day] for day in days) == 3 * y50[course] + 2 * y80[course],
                name=f"TotalDays_{course}"
            )

            for i in range(len(days)):
                for j in range(i + 1, len(days)):
                    day_i = days[i]
                    day_j = days[j]
                    day_num_i = day_numbers[day_i]
                    day_num_j = day_numbers[day_j]
                    # If days are consecutive
                    if abs(day_num_i - day_num_j) == 1:
                        # Prevent both consecutive days from being chosen simultaneously
                        model.addConstr(
                            d[course][day_i] + d[course][day_j] <= 1,
                            name=f"NoConsecutiveDays_{course}_{day_i}_{day_j}"
                        )

            # Ensuring consistent start times for non-lab courses:
            h_50 = model.addVars(distinct_start_times_50, vtype=GRB.BINARY, name=f"h50_{course}")
            h_80 = model.addVars(distinct_start_times_80, vtype=GRB.BINARY, name=f"h80_{course}")

            model.addConstr(
                gp.quicksum(h_50[t] for t in distinct_start_times_50) == y50[course],
                name=f"OneStartTime50_{course}"
            )
            model.addConstr(
                gp.quicksum(h_80[t] for t in distinct_start_times_80) == y80[course],
                name=f"OneStartTime80_{course}"
            )

            # Link x[c,s] to h_...:
            for t in distinct_start_times_50:
                for ss in slots_for_start_50[t]:
                    model.addConstr(
                        x[course, ss] <= h_50[t],
                        name=f"LinkXH50_{course}_{ss}"
                    )
            for t in distinct_start_times_80:
                for ss in slots_for_start_80[t]:
                    model.addConstr(
                        x[course, ss] <= h_80[t],
                        name=f"LinkXH80_{course}_{ss}"
                    )

# consider for discussion
# Add non-overlapping constraints for courses in the same year from flowchart
for year, course_list in flowchart.items():
    # Filter to offered courses
    filtered_list = []
    for c in course_list:
        if c in lab_courses_dict:
            for i in range(1, lab_courses_dict[c]['sections'] + 1):
                filtered_list.append(f"{c}_sec{i}")
        elif c in courses:
            filtered_list.append(c)

    if len(filtered_list) > 1:
        # Original line
        for i in range(len(filtered_list)):
            for j in range(i+1, len(filtered_list)):
                c1 = filtered_list[i]
                c2 = filtered_list[j]
                if c1 in courses and c2 in courses:
                    for s in slots:
                        model.addConstr(
                            x[c1, s] + x[c2, s] <= 1,
                            name=f"NoCoreqSameSlot_{c1}_{c2}_{s}"
                        )
                    for s1 in slots:
                        if s1 in overlaps:
                            for s2 in overlaps[s1]:
                                model.addConstr(
                                    x[c1, s1] + x[c2, s2] <= 1,
                                    name=f"FlowchartConflict_{year}_{c1}_{c2}_{s1}_{s2}"
                                )

        # Discussion line
        discussion_lines = create_discussion_line(filtered_list)
        for discussion_line in discussion_lines:
            if len(discussion_line) > 1:
                for i in range(len(discussion_line)):
                    for j in range(i+1, len(discussion_line)):
                        c1 = discussion_line[i]
                        c2 = discussion_line[j]
                        if c1 in courses and c2 in courses:
                            for s in slots:
                                model.addConstr(
                                    x[c1, s] + x[c2, s] <= 1,
                                    name=f"NoCoreqSameSlot_{c1}_{c2}_{s}"
                                )
                            for s1 in slots:
                                if s1 in overlaps:
                                    for s2 in overlaps[s1]:
                                        model.addConstr(
                                            x[c1, s1] + x[c2, s2] <= 1,
                                            name=f"FlowchartConflict_{year}_{c1}_{c2}_{s1}_{s2}"
                                        )
# consider for discussion
# Corequisite constraints
for index, row in df_corequisites.iterrows():
    course1 = row['1']
    course2 = row['2']
    pair = []
    if course1 in courses:
        pair.append(course1)
    if course2 in courses:
        pair.append(course2)

    if len(pair) == 2:
        # Original pair line

        # 1) No same-slot
        for s in slots:
            model.addConstr(
                x[pair[0], s] + x[pair[1], s] <= 1,
                name=f"NoCoreqSameSlot_{pair[0]}_{pair[1]}_{s}"
            )
        for s1 in slots:
            if s1 in overlaps:
                for s2 in overlaps[s1]:
                    model.addConstr(
                        x[pair[0], s1] + x[pair[1], s2] <= 1,
                        name=f"CoreqConflict_{pair[0]}_{pair[1]}_{s1}_{s2}"
                    )

        # Discussion line
        disc_pairs = create_discussion_line(pair)
        for disc_pair in disc_pairs:
            if len(disc_pair) == 2:
                for s in slots:
                    model.addConstr(
                        x[disc_pair[0], s] + x[disc_pair[1], s] <= 1,
                        name=f"NoCoreqSameSlot_{disc_pair[0]}_{disc_pair[1]}_{s}"
                    )

                for s1 in slots:
                    if s1 in overlaps:
                        for s2 in overlaps[s1]:
                            model.addConstr(
                                x[disc_pair[0], s1] + x[disc_pair[1], s2] <= 1,
                                name=f"CoreqConflict_{disc_pair[0]}_{disc_pair[1]}_{s1}_{s2}"
                            )


# Penalties for early slots and Friday late-slots
for c in courses:
    for s in penalty_slots:
        model.addConstr(w_penalty[c, s] == x[c, s], name=f"EarlyorFridayLateSlot_{c}_{s}")


# Overlap penalties for 3000-level courses
for i in range(len(courses_3000)):
    for j in range(i+1, len(courses_3000)):
        c1 = courses_3000[i]
        c2 = courses_3000[j]
        for s1 in slots:
            if s1 in overlaps:
                for s2 in overlaps[s1]:
                    overlap_var = model.addVar(vtype=GRB.BINARY, name=f"Overlap3000_{c1}_{c2}_{s1}_{s2}")
                    model.addConstr(overlap_var >= x[c1, s1] + x[c2, s2] - 1)
                    # Store in w_overlap_3000
                    w_overlap_3000[(c1, c2, s1, s2)] = overlap_var

# Overlap penalties for design vs 3000-level or capstone
for c_design in design_courses:
    if c_design in courses:
        # With 3000-level courses
        for c_3000 in courses_3000:
            if c_3000 in courses:
                for s1 in slots:
                    if s1 in overlaps:
                        for s2 in overlaps[s1]:
                            overlap_var = model.addVar(vtype=GRB.BINARY, name=f"OverlapDesign_{c_design}_{c_3000}_{s1}_{s2}")
                            model.addConstr(overlap_var >= x[c_design, s1] + x[c_3000, s2] - 1)
                            w_overlap_design[(c_design, c_3000, s1, s2)] = overlap_var

        # With capstone courses
        for c_cap in capstone_courses:
            if c_cap in courses:
                for s1 in slots:
                    if s1 in overlaps:
                        for s2 in overlaps[s1]:
                            overlap_var = model.addVar(vtype=GRB.BINARY, name=f"OverlapDesignCapstone_{c_design}_{c_cap}_{s1}_{s2}")
                            model.addConstr(overlap_var >= x[c_design, s1] + x[c_cap, s2] - 1)
                            # Add to w_overlap_design as well since design vs capstone also undesired
                            w_overlap_design[(c_design, c_cap, s1, s2)] = overlap_var

# Overlap penalties for tech electives vs capstone
for c_tech in tech_electives:
    if c_tech in courses:
        for c_cap in capstone_courses:
            if c_cap in courses:
                for s1 in slots:
                    if s1 in overlaps:
                        for s2 in overlaps[s1]:
                            overlap_var = model.addVar(vtype=GRB.BINARY, name=f"OverlapTech_{c_tech}_{c_cap}_{s1}_{s2}")
                            model.addConstr(overlap_var >= x[c_tech, s1] + x[c_cap, s2] - 1)
                            w_overlap_tech[(c_tech, c_cap, s1, s2)] = overlap_var

# Add constraints involving discussion
for dcourse, req in discussion_courses.items():
    d_disc = model.addVars(days, vtype=GRB.BINARY, name=f"d_{dcourse}")

    if req['type'] == 'mixed':
        y50_disc = model.addVar(vtype=GRB.BINARY, name=f"y50_disc_{dcourse.replace(' ','_')}")
        y80_disc = model.addVar(vtype=GRB.BINARY, name=f"y80_disc_{dcourse.replace(' ','_')}")
        model.addConstr(y50_disc + y80_disc == 1, name=f"DiscSlotLengthChoice_{dcourse}")

        # Find sequences for mixed discussion course (like a lab)
        sequence_vars = []
        sequence_info = []
        for day in days:
            # For 50-min sequences
            slots_50_day = [ss for ss in slots if slot_days[ss] == day and slot_lengths[ss] == 50]
            slots_50_day_sorted = sorted(slots_50_day, key=lambda ss: slot_start_mins[ss])
            for idx in range(len(slots_50_day_sorted)-req['50_min_needed']+1):
                seq_slots = slots_50_day_sorted[idx: idx+req['50_min_needed']]
                times = [slot_start_mins[s] for s in seq_slots]
                if all(times[i+1]-times[i] == 55 for i in range(len(times)-1)):
                    seq_var = model.addVar(vtype=GRB.BINARY, name=f"seq_{dcourse}_{day}_50_{idx}")
                    sequence_vars.append(seq_var)
                    sequence_info.append({'day': day, 'length': 50, 'slots': seq_slots, 'var': seq_var})

            # For 80-min sequences
            slots_80_day = [ss for ss in slots if slot_days[ss] == day and slot_lengths[ss] == 80]
            slots_80_day_sorted = sorted(slots_80_day, key=lambda ss: slot_start_mins[ss])
            for idx in range(len(slots_80_day_sorted)-req['80_min_needed']+1):
                seq_slots = slots_80_day_sorted[idx: idx+req['80_min_needed']]
                times = [slot_start_mins[s] for s in seq_slots]
                if all(times[i+1]-times[i] == 95 for i in range(len(times)-1)):
                    seq_var = model.addVar(vtype=GRB.BINARY, name=f"seq_{dcourse}_{day}_80_{idx}")
                    sequence_vars.append(seq_var)
                    sequence_info.append({'day': day, 'length': 80, 'slots': seq_slots, 'var': seq_var})

        model.addConstr(gp.quicksum(sequence_vars) == 1, name=f"DiscSeq_{dcourse}")

        for seq in sequence_info:
            seq_var = seq['var']
            seq_slots = seq['slots']
            length = seq['length']
            if length == 50:
                model.addConstr(seq_var <= y50_disc, name=f"SeqVar50_{dcourse}")
            else:
                model.addConstr(seq_var <= y80_disc, name=f"SeqVar80_{dcourse}")
            for s in seq_slots:
                model.addConstr(x[dcourse, s] == seq_var, name=f"DiscSeqAssign_{dcourse}_{s}")

        model.addConstr(
            gp.quicksum(x[dcourse, s] for s in slots) == 3*y50_disc + 2*y80_disc,
            name=f"DiscTotalSlots_{dcourse}"
        )

        for day in days:
            seqs_day = [seq['var'] for seq in sequence_info if seq['day']==day]
            model.addConstr(
                d_disc[day] == gp.quicksum(seqs_day),
                name=f"DiscDay_{dcourse}_{day}"
            )
        model.addConstr(
            gp.quicksum(d_disc[day] for day in days) == 1,
            name=f"DiscSingleDay_{dcourse}"
        )

    elif req['type'] == '50_only':
        y50_disc = model.addVar(vtype=GRB.BINARY, name=f"y50_disc_{dcourse.replace(' ','_')}")
        model.addConstr(
            gp.quicksum(x[dcourse, s] for s in slots if slot_lengths[s]==50) == y50_disc,
            name=f"{dcourse}_One50Slot"
        )
        # If main course offered => discussion offered, so y50_disc=1
        model.addConstr(y50_disc == 1, name=f"DiscOffered_{dcourse}")

        for day in days:
            model.addConstr(
                d_disc[day] == gp.quicksum(x[dcourse, s] for s in slots if slot_days[s]==day),
                name=f"DiscDay_{dcourse}_{day}"
            )
        model.addConstr(
            gp.quicksum(d_disc[day] for day in days) == 1,
            name=f"DiscSingleDay_{dcourse}"
        )

    elif req['type'] == '80_only':
        y80_disc = model.addVar(vtype=GRB.BINARY, name=f"y80_disc_{dcourse.replace(' ','_')}")
        model.addConstr(
            gp.quicksum(x[dcourse, s] for s in slots if slot_lengths[s]==80) == y80_disc,
            name=f"{dcourse}_One80Slot"
        )
        model.addConstr(y80_disc == 1, name=f"DiscOffered_{dcourse}")

        for day in days:
            model.addConstr(
                d_disc[day] == gp.quicksum(x[dcourse, s] for s in slots if slot_days[s]==day),
                name=f"DiscDay_{dcourse}_{day}"
            )
        model.addConstr(
            gp.quicksum(d_disc[day] for day in days) == 1,
            name=f"DiscSingleDay_{dcourse}"
        )


# Add constraints to prevent overlap between main course and its discussion course
for dcourse in discussion_courses.keys():
    main_course = discussion_to_main[dcourse]

    # First, prevent them from choosing the exact same slot:
    for s in slots:
        # If both main course and discussion course choose the same slot s:
        model.addConstr(
            x[main_course, s] + x[dcourse, s] <= 1,
            name=f"NoMainDiscSameSlot_{main_course}_{dcourse}_{s}"
        )

    # Now, also prevent them from choosing different overlapping slots:
    for s1 in slots:
        if s1 in overlaps:
            for s2 in overlaps[s1]:
                # If s1 and s2 overlap in time, prevent main_course and dcourse from using them simultaneously
                model.addConstr(
                    x[main_course, s1] + x[dcourse, s2] <= 1,
                    name=f"NoMainDiscOverlap1_{main_course}_{dcourse}_{s1}_{s2}"
                )
                model.addConstr(
                    x[dcourse, s1] + x[main_course, s2] <= 1,
                    name=f"NoMainDiscOverlap2_{main_course}_{dcourse}_{s1}_{s2}"
                )

# Imposing No-Overlap Constraints among discussion sections
# For discussion courses:
for main_c, (disc_c, disc_type, disc_num_secs) in required_discussions.items():
    if main_c in courses_offered:
        # We have discussion_courses_map[disc_c] storing info
        sec_names = discussion_courses_map[disc_c]['sections']
        # Imposing no overlap among these discussion sections
        for i in range(len(sec_names)):
            for j in range(i+1, len(sec_names)):
                sec1 = sec_names[i]
                sec2 = sec_names[j]
                # No same slot
                for s in slots:
                    model.addConstr(
                        x[sec1, s] + x[sec2, s] <= 1,
                        name=f"NoSameSlot_{sec1}_{sec2}_{s}"
                    )
                # No overlapping distinct slots
                for s1 in slots:
                    if s1 in overlaps:
                        for s2 in overlaps[s1]:
                            model.addConstr(
                                x[sec1, s1] + x[sec2, s2] <= 1,
                                name=f"NoOverlap_{sec1}_{sec2}_{s1}_{s2}"
                            )

# Objective weights
penalty_slot = -1
penalty_overlap_3000 = -1
penalty_overlap_design = -1
penaly_overlap_tech = -1

model.setObjective(
    gp.quicksum(z[c] for c in courses) +
    gp.quicksum(u_var for u_var in u.values()) +
    penalty_slot * gp.quicksum(w_penalty[c, s] for c in courses for s in penalty_slots) +
    penalty_overlap_3000 * gp.quicksum(w_overlap_3000.values()) +
    penalty_overlap_design * gp.quicksum(w_overlap_design.values()) +
    penaly_overlap_tech * gp.quicksum(w_overlap_tech.values()),
    GRB.MAXIMIZE
)

model.optimize()

# After model.optimize() and ensuring multiple solutions are found
if model.SolCount == 0:
    print("No solutions found.")
else:
    numSolutions = min(model.SolCount, 3)

    with pd.ExcelWriter(plan_semester + "_schedule.xlsx", engine="xlsxwriter") as writer:
        model.Params.SolutionNumber = 0

        schedule = []
        for course in courses:
            assigned_slots = [s for s in slots if x[course, s].Xn > 0.5]
            for s in assigned_slots:
                day = df_slots.loc[df_slots['SlotID'] == s, 'Day'].values[0]
                start_time = df_slots.loc[df_slots['SlotID'] == s, 'Start'].values[0]
                end_time = df_slots.loc[df_slots['SlotID'] == s, 'End'].values[0]
                schedule.append({
                    'Course': course,
                    'Day': day,
                    'Start': start_time,
                    'End': end_time
                })

        df_schedule = pd.DataFrame(schedule)

        # Combine slots on the same day
        df_day_combined = (df_schedule
                            .groupby(['Course','Day'], as_index=False)
                            .agg({'Start':'min','End':'max'}))

        # Combine multiple days for each course
        df_final = (df_day_combined
                    .groupby('Course', as_index=False)[['Day','Start','End']]
                    .apply(lambda g: pd.Series({
                        'Day': ",".join([day_short[d] for d in sorted(g['Day'], key=lambda d: day_numbers[d])]),
                        'Start': g['Start'].min(),
                        'End': g['End'].max()
                    }))
                    .reset_index(drop=True))

        df_final = df_final[['Course','Day','Start','End']]

        sheet_name = f"Optimal_Solution"
        df_final.to_excel(writer, sheet_name=sheet_name, index=False)

        # Access the workbook and worksheet objects
        workbook  = writer.book
        worksheet = writer.sheets[sheet_name]

        # Adjust column widths
        for idx, col in enumerate(df_final.columns):
            max_len = max(
                df_final[col].astype(str).map(len).max(),
                len(col)
            )
            worksheet.set_column(idx, idx, max_len + 2)