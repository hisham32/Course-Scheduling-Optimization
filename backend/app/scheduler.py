import os
import re
import io
import pandas as pd
import numpy as np
from datetime import datetime
from typing import Dict, Any, Optional, Tuple

import gurobipy as gp
from gurobipy import GRB


def _time_to_minutes(t: str) -> int:
    h, m, s = map(int, t.split(":"))
    return h * 60 + m


def _extract_first_duration(meeting_pattern: str) -> Optional[int]:
    if pd.isna(meeting_pattern) or "Does Not Meet" in str(meeting_pattern):
        return None
    try:
        time_match = re.search(r'(\d{1,2}(:\d{2})?[ap]m)-(\d{1,2}(:\d{2})?[ap]m)', str(meeting_pattern))
        if time_match:
            start_time, _, end_time, _ = time_match.groups()
            if ":" not in start_time:
                start_time = start_time.replace("am", ":00am").replace("pm", ":00pm")
            if ":" not in end_time:
                end_time = end_time.replace("am", ":00am").replace("pm", ":00pm")
            start_time = pd.to_datetime(start_time.strip(), format="%I:%M%p")
            end_time = pd.to_datetime(end_time.strip(), format="%I:%M%p")
            return int((end_time - start_time).seconds // 60)
    except Exception:
        return None
    return None


def _read_env_wls_params() -> Dict[str, Any]:
    wls_id = os.getenv("GRB_WLSACCESSID")
    wls_secret = os.getenv("GRB_WLSSECRET")
    license_id = os.getenv("GRB_LICENSEID")
    if wls_id and wls_secret and license_id:
        try:
            license_id = int(license_id)
        except ValueError:
            pass
        return {"WLSACCESSID": wls_id, "WLSSECRET": wls_secret, "LICENSEID": license_id}
    return {}


def extract_latest_preferences(offerings_path: str) -> pd.DataFrame:
    prev_offerings = pd.ExcelFile(offerings_path)
    semester_order = prev_offerings.sheet_names

    all_data = []
    for sheet in semester_order:
        df = prev_offerings.parse(sheet)
        if df.empty:
            continue
        df["Semester"] = sheet
        df["Course"] = df["Subject Code"] + " " + df["Catalog Number"].astype(str)
        all_data.append(df[["Course", "Meeting Pattern", "Semester"]])

    combined_data = pd.concat(all_data, ignore_index=True)
    combined_data["Semester"] = pd.Categorical(combined_data["Semester"], categories=semester_order, ordered=True)
    combined_data = combined_data.sort_values(by=["Course", "Semester"], ascending=[True, False])

    combined_data = combined_data.sort_values(
        by=["Course", "Meeting Pattern"],
        key=lambda x: x != "Does Not Meet"
    ).drop_duplicates(subset=["Course"], keep="last")

    latest_data = combined_data.groupby("Course").first().reset_index()
    latest_data["Duration"] = latest_data["Meeting Pattern"].apply(_extract_first_duration)
    latest_data["50_min"] = (latest_data["Duration"] == 50).astype(int)
    latest_data["80_min"] = (latest_data["Duration"] == 80).astype(int)
    return latest_data[["Course", "50_min", "80_min"]].copy()


def run_optimization(plan_semester: str, offerings_path: str, input_data_path: str) -> Tuple[pd.DataFrame, str]:
    # Create Gurobi env (WLS if available)
    wls_params = _read_env_wls_params()
    env = gp.Env(params=wls_params) if wls_params else gp.Env(empty=True)

    # Load inputs
    catalog_df      = pd.read_excel(input_data_path, sheet_name="catalog")
    flowchart_df    = pd.read_excel(input_data_path, sheet_name="flowchart")
    timeslots50_df  = pd.read_excel(input_data_path, sheet_name="timeslots50")
    timeslots80_df  = pd.read_excel(input_data_path, sheet_name="timeslots80")
    df_prerequisites= pd.read_excel(input_data_path, sheet_name="prerequisites")
    df_corequisites = pd.read_excel(input_data_path, sheet_name="corequisites")
    df_capstone     = pd.read_excel(input_data_path, sheet_name="capstone")
    df_design       = pd.read_excel(input_data_path, sheet_name="design")
    df_tech         = pd.read_excel(input_data_path, sheet_name="tech")
    df_discussion   = pd.read_excel(input_data_path, sheet_name="discussion", index_col=0)
    df_lab          = pd.read_excel(input_data_path, sheet_name="lab", index_col=0)

    catalog_df["course"] = catalog_df["Subject code"] + " " + catalog_df["Catalog Number"].astype(str)
    courses_offered = []
    for _, row in catalog_df.iterrows():
        course = row["course"]
        if course == "CVEEN 5920":
            continue
        if row["Semesters Typically Offered"] == plan_semester or row["Semesters Typically Offered"] == "All":
            courses_offered.append(course)

    df_preference = extract_latest_preferences(offerings_path)

    slots_50_min = timeslots50_df.to_dict(orient="list")
    slots_80_min = timeslots80_df.to_dict(orient="list")
    flowchart = flowchart_df.to_dict(orient="list")

    # Build slots
    slots_50 = {"Day": [], "Start": [], "End": []}
    for day in ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday"]:
        start_key = f"{day}_start"
        end_key = f"{day}_end"
        if start_key in slots_50_min:
            start_times = [t for t in slots_50_min[start_key] if pd.notna(t)]
            end_times = [t for t in slots_50_min[end_key] if pd.notna(t)]
            slots_50["Day"].extend([day] * len(start_times))
            slots_50["Start"].extend(start_times)
            slots_50["End"].extend(end_times)
    df_slots_50 = pd.DataFrame(slots_50)

    slots_80 = {"Day": [], "Start": [], "End": []}
    for day in ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday"]:
        start_key = f"{day}_start"
        end_key = f"{day}_end"
        if start_key in slots_80_min:
            start_times = [t for t in slots_80_min[start_key] if pd.notna(t)]
            end_times = [t for t in slots_80_min[end_key] if pd.notna(t)]
            slots_80["Day"].extend([day] * len(start_times))
            slots_80["Start"].extend(start_times)
            slots_80["End"].extend(end_times)
    df_slots_80 = pd.DataFrame(slots_80)

    df_slots_50["SlotLength"] = 50
    df_slots_80["SlotLength"] = 80
    df_slots = pd.concat([df_slots_50, df_slots_80], ignore_index=True)
    df_slots["SlotID"] = df_slots.index
    df_slots["StartMin"] = df_slots["Start"].apply(_time_to_minutes)
    df_slots["EndMin"] = df_slots["End"].apply(_time_to_minutes)

    # Model
    model = gp.Model("Course Scheduling", env=env)

    # Configure Gurobi log file (truncate first for fresh run)
    LOG_PATH = "/tmp/gurobi.log"
    try:
        with open(LOG_PATH, "w", encoding="utf-8") as _f:
            pass
    except Exception:
        pass
    model.setParam("LogFile", LOG_PATH)

    slots = df_slots["SlotID"].tolist()
    days = df_slots["Day"].unique().tolist()

    day_numbers = {"Monday": 1, "Tuesday": 2, "Wednesday": 3, "Thursday": 4, "Friday": 5}
    day_short = {"Monday": "Mo", "Tuesday": "Tu", "Wednesday": "We", "Thursday": "Th", "Friday": "Fr"}

    slot_lengths = df_slots.set_index("SlotID")["SlotLength"].to_dict()
    slot_days = df_slots.set_index("SlotID")["Day"].to_dict()
    slot_start_times = df_slots.set_index("SlotID")["Start"].to_dict()
    slot_start_mins = df_slots.set_index("SlotID")["StartMin"].to_dict()
    slot_end_mins = df_slots.set_index("SlotID")["EndMin"].to_dict()

    capstone_courses = df_capstone["course"].tolist()
    design_courses = df_design["course"].tolist()
    tech_electives = df_tech["course"].tolist()

    lab_courses_dict = df_lab.to_dict(orient="dict")

    courses = list(set([c for c in courses_offered if (c not in lab_courses_dict) and (1000 <= int(c.split()[1]) < 5999)]))
    courses_3000 = [c for c in courses_offered if (c not in lab_courses_dict) and (3000 <= int(c.split()[1]) < 4000)]

    # Lab sections
    lab_sections_map = {}
    lab_courses = []
    for main_lab, num_secs in lab_courses_dict.items():
        section_names = []
        for i in range(1, num_secs["sections"] + 1):
            sec_name = f"{main_lab}_sec{i}"
            section_names.append(sec_name)
            courses.append(sec_name)
            lab_courses.append(sec_name)
            courses_3000.append(sec_name)
        lab_sections_map[main_lab] = section_names

    # Discussions
    required_discussions = df_discussion.to_dict(orient="list")

    def get_discussion_course(course):
        for main_c, (disc_c, ctype, sections) in required_discussions.items():
            if main_c == course:
                return [f"{disc_c}_sec{i}" for i in range(1, sections + 1)]
        return None

    def create_discussion_line(course_list):
        discussion_lines = []
        replaced = False
        for i, c in enumerate(course_list):
            dcourse = get_discussion_course(c)
            if not replaced and dcourse is not None and all(sec in courses for sec in dcourse):
                for disc_sec in dcourse:
                    new_line = course_list.copy()
                    new_line[i] = disc_sec
                    discussion_lines.append(new_line)
                replaced = True
        if not replaced:
            discussion_lines.append(course_list)
        return discussion_lines

    discussion_to_main = {}
    for main_c, (disc_c, ctype, sections) in required_discussions.items():
        for i in range(1, sections + 1):
            discussion_to_main[f"{disc_c}_sec{i}"] = main_c

    discussion_courses = {}
    for main_course, (disc_course, ctype, sections) in required_discussions.items():
        if main_course in courses_offered:
            for i in range(1, sections + 1):
                disc_course_sec = f"{disc_course}_sec{i}"
                discussion_courses[disc_course_sec] = {
                    "type": ctype,
                    "50_min_needed": 3 if ctype == "mixed" else (1 if ctype == "50_only" else 0),
                    "80_min_needed": 2 if ctype == "mixed" else (1 if ctype == "80_only" else 0),
                }

    discussion_courses_map = {}
    for main_c, (disc_c, disc_type, disc_num_secs) in required_discussions.items():
        if main_c in courses_offered:
            sec_names = []
            for i in range(1, disc_num_secs + 1):
                sec_name = f"{disc_c}_sec{i}"
                sec_names.append(sec_name)
                courses.append(sec_name)
            discussion_courses_map[disc_c] = {"type": disc_type, "num_secs": disc_num_secs, "sections": sec_names}

    # Variables
    x = model.addVars(courses, slots, vtype=GRB.BINARY, name="x")
    y50 = model.addVars(courses, vtype=GRB.BINARY, name="y50")
    y80 = model.addVars(courses, vtype=GRB.BINARY, name="y80")
    z = model.addVars(courses, vtype=GRB.BINARY, name="z")

    distinct_start_times_50 = sorted(set(df_slots.loc[df_slots["SlotLength"] == 50, "StartMin"]))
    distinct_start_times_80 = sorted(set(df_slots.loc[df_slots["SlotLength"] == 80, "StartMin"]))

    slots_for_start_50 = {t: [] for t in distinct_start_times_50}
    slots_for_start_80 = {t: [] for t in distinct_start_times_80}
    for s in slots:
        if df_slots.loc[s, "SlotLength"] == 50:
            slots_for_start_50[df_slots.loc[s, "StartMin"]].append(s)
        else:
            slots_for_start_80[df_slots.loc[s, "StartMin"]].append(s)

    # Overlaps
    overlaps = {}
    for s1 in slots:
        for s2 in slots:
            if s1 >= s2:
                continue
            if df_slots.loc[s1, "Day"] != df_slots.loc[s2, "Day"]:
                continue
            start1, end1 = df_slots.loc[s1, "StartMin"], df_slots.loc[s1, "EndMin"]
            start2, end2 = df_slots.loc[s2, "StartMin"], df_slots.loc[s2, "EndMin"]
            if max(start1, start2) < min(end1, end2):
                overlaps.setdefault(s1, set()).add(s2)
                overlaps.setdefault(s2, set()).add(s1)

    # Prereqs + discussion variants
    u = {}
    for _, row in df_prerequisites.iterrows():
        prereq_courses = [row[str(i)] for i in range(1, len(row)) if pd.notna(row.get(str(i)))]
        prereq_courses = [c for c in prereq_courses if c in courses]
        if len(prereq_courses) <= 1:
            continue
        for s in slots:
            u_vars = [x[c, s] for c in prereq_courses]
            u_var = model.addVar(vtype=GRB.BINARY, name=f"Prereq_{'_'.join(prereq_courses)}_{s}")
            model.addGenConstrMin(u_var, u_vars)
            u[(tuple(prereq_courses), s)] = u_var

        discussion_lines = create_discussion_line(prereq_courses)
        for discussion_line in discussion_lines:
            if len(discussion_line) > 1:
                for s in slots:
                    u_vars = [x[c, s] for c in discussion_line]
                    u_var = model.addVar(vtype=GRB.BINARY, name=f"Prereq_{'_'.join(discussion_line)}_{s}")
                    model.addGenConstrMin(u_var, u_vars)
                    u[(tuple(discussion_line), s)] = u_var

    # Early/Friday-late penalties
    early_slots = [s for s in slots if df_slots.loc[s, "Start"] == "07:30:00"]
    friday_late_slots = [s for s in slots if df_slots.loc[s, "Day"] == "Friday" and df_slots.loc[s, "StartMin"] >= 14 * 60]
    penalty_slots = early_slots + friday_late_slots
    w_penalty = model.addVars(courses, penalty_slots, vtype=GRB.BINARY, name="wFriLate")

    w_overlap_3000 = {}
    w_overlap_design = {}
    w_overlap_tech = {}

    # Per-course constraints (same logic as before)
    for course in courses:
        if course in discussion_courses.keys():
            continue

        if (course in [*discussion_courses.keys()]) or (course not in df_preference["Course"].tolist()):
            pref_50 = 1
            pref_80 = 1
        else:
            pref_50 = int(df_preference.loc[df_preference["Course"] == course, "50_min"].values[0])
            pref_80 = int(df_preference.loc[df_preference["Course"] == course, "80_min"].values[0])

        if course in [*lab_courses]:
            model.addConstr(z[course] == 0)
            model.addConstr(y50[course] + y80[course] == 1)
        elif course in df_capstone["course"].tolist():
            model.addConstr(z[course] == 0)
            model.addConstr(y50[course] == 0)
            model.addConstr(y80[course] == 1)
        else:
            model.addConstr(y50[course] + y80[course] == 1)
            if pref_50 == 1 and pref_80 == 0:
                model.addConstr(z[course] <= y50[course])
            elif pref_80 == 1 and pref_50 == 0:
                model.addConstr(z[course] <= y80[course])
            else:
                model.addConstr(z[course] <= y50[course] + y80[course])

        if course in lab_courses:
            n_slots_50 = 3
            n_slots_80 = 2
            sequence_vars = []
            sequence_info = []
            for day in days:
                # 50
                slots_50_day = [ss for ss in slots if slot_days[ss] == day and slot_lengths[ss] == 50]
                slots_50_day_sorted = sorted(slots_50_day, key=lambda ss: slot_start_mins[ss])
                for idx in range(len(slots_50_day_sorted) - n_slots_50 + 1):
                    seq_slots = slots_50_day_sorted[idx: idx + n_slots_50]
                    times = [slot_start_mins[ss] for ss in seq_slots]
                    if all(times[i + 1] - times[i] == 55 for i in range(len(times) - 1)):
                        seq_var = model.addVar(vtype=GRB.BINARY, name=f"seq_{course}_{day}_50_{idx}")
                        sequence_vars.append(seq_var)
                        sequence_info.append({"day": day, "length": 50, "slots": seq_slots, "var": seq_var})
                # 80
                slots_80_day = [ss for ss in slots if slot_days[ss] == day and slot_lengths[ss] == 80]
                slots_80_day_sorted = sorted(slots_80_day, key=lambda ss: slot_start_mins[ss])
                for idx in range(len(slots_80_day_sorted) - n_slots_80 + 1):
                    seq_slots = slots_80_day_sorted[idx: idx + n_slots_80]
                    times = [slot_start_mins[ss] for ss in seq_slots]
                    if all(times[i + 1] - times[i] == 95 for i in range(len(times) - 1)):
                        seq_var = model.addVar(vtype=GRB.BINARY, name=f"seq_{course}_{day}_80_{idx}")
                        sequence_vars.append(seq_var)
                        sequence_info.append({"day": day, "length": 80, "slots": seq_slots, "var": seq_var})

            model.addConstr(gp.quicksum(sequence_vars) == 1)
            for seq in sequence_info:
                for ss in seq["slots"]:
                    model.addConstr(x[course, ss] == seq["var"])
            model.addConstr(
                gp.quicksum(x[course, ss] for ss in slots) ==
                gp.quicksum(len(seq["slots"]) * seq["var"] for seq in sequence_info)
            )
            d = model.addVars(days, vtype=GRB.BINARY, name=f"d_{course}")
            for day in days:
                model.addConstr(d[day] == gp.quicksum(seq["var"] for seq in sequence_info if seq["day"] == day))
            model.addConstr(gp.quicksum(d[day] for day in days) == 1)

        else:
            if course == "CVEEN 2000":
                model.addConstr(y50[course] == 1)
                model.addConstr(y80[course] == 0)
                model.addConstr(gp.quicksum(x[course, s] for s in slots if slot_lengths[s] == 50) == 1)
            else:
                model.addConstr(gp.quicksum(x[course, ss] for ss in slots if slot_lengths[ss] == 50) == 3 * y50[course])
                model.addConstr(gp.quicksum(x[course, ss] for ss in slots if slot_lengths[ss] == 80) == 2 * y80[course])
                for day in days:
                    model.addConstr(gp.quicksum(x[course, ss] for ss in slots if slot_days[ss] == day) <= 1)

                d_course = {day: model.addVar(vtype=GRB.BINARY, name=f"d_{course}_{day}") for day in days}
                for day in days:
                    slots_on_day = [ss for ss in slots if slot_days[ss] == day]
                    model.addConstr(d_course[day] == gp.quicksum(x[course, ss] for ss in slots_on_day))
                model.addConstr(gp.quicksum(d_course[day] for day in days) == 3 * y50[course] + 2 * y80[course])

                for i in range(len(days)):
                    for j in range(i + 1, len(days)):
                        day_i = days[i]
                        day_j = days[j]
                        if abs(day_numbers[day_i] - day_numbers[day_j]) == 1:
                            model.addConstr(d_course[day_i] + d_course[day_j] <= 1)

                h_50 = {t: model.addVar(vtype=GRB.BINARY, name=f"h50_{course}_{t}") for t in distinct_start_times_50}
                h_80 = {t: model.addVar(vtype=GRB.BINARY, name=f"h80_{course}_{t}") for t in distinct_start_times_80}
                model.addConstr(gp.quicksum(h_50[t] for t in distinct_start_times_50) == y50[course])
                model.addConstr(gp.quicksum(h_80[t] for t in distinct_start_times_80) == y80[course])
                for t in distinct_start_times_50:
                    for ss in slots_for_start_50[t]:
                        model.addConstr(x[course, ss] <= h_50[t])
                for t in distinct_start_times_80:
                    for ss in slots_for_start_80[t]:
                        model.addConstr(x[course, ss] <= h_80[t])

    # Flowchart conflicts (and discussion variants)
    for year, course_list in flowchart.items():
        filtered_list = []
        for c in course_list:
            if c in lab_courses_dict:
                for i in range(1, lab_courses_dict[c]["sections"] + 1):
                    filtered_list.append(f"{c}_sec{i}")
            elif c in courses:
                filtered_list.append(c)
        if len(filtered_list) > 1:
            for i in range(len(filtered_list)):
                for j in range(i + 1, len(filtered_list)):
                    c1, c2 = filtered_list[i], filtered_list[j]
                    if c1 in courses and c2 in courses:
                        for s in slots:
                            model.addConstr(x[c1, s] + x[c2, s] <= 1)
                        for s1 in slots:
                            if s1 in overlaps:
                                for s2 in overlaps[s1]:
                                    model.addConstr(x[c1, s1] + x[c2, s2] <= 1)
            for discussion_line in create_discussion_line(filtered_list):
                if len(discussion_line) > 1:
                    for i in range(len(discussion_line)):
                        for j in range(i + 1, len(discussion_line)):
                            c1, c2 = discussion_line[i], discussion_line[j]
                            if c1 in courses and c2 in courses:
                                for s in slots:
                                    model.addConstr(x[c1, s] + x[c2, s] <= 1)
                                for s1 in slots:
                                    if s1 in overlaps:
                                        for s2 in overlaps[s1]:
                                            model.addConstr(x[c1, s1] + x[c2, s2] <= 1)

    # Coreqs (+ discussion variants)
    for _, row in df_corequisites.iterrows():
        course1 = row["1"]
        course2 = row["2"]
        pair = []
        if course1 in courses:
            pair.append(course1)
        if course2 in courses:
            pair.append(course2)
        if len(pair) == 2:
            for s in slots:
                model.addConstr(x[pair[0], s] + x[pair[1], s] <= 1)
            for s1 in slots:
                if s1 in overlaps:
                    for s2 in overlaps[s1]:
                        model.addConstr(x[pair[0], s1] + x[pair[1], s2] <= 1)
            for disc_pair in create_discussion_line(pair):
                if len(disc_pair) == 2:
                    for s in slots:
                        model.addConstr(x[disc_pair[0], s] + x[disc_pair[1], s] <= 1)
                    for s1 in slots:
                        if s1 in overlaps:
                            for s2 in overlaps[s1]:
                                model.addConstr(x[disc_pair[0], s1] + x[disc_pair[1], s2] <= 1)

    # Penalties
    for c in courses:
        for s in penalty_slots:
            model.addConstr(w_penalty[c, s] == x[c, s])

    for i in range(len(courses_3000)):
        for j in range(i + 1, len(courses_3000)):
            c1, c2 = courses_3000[i], courses_3000[j]
            for s1 in slots:
                if s1 in overlaps:
                    for s2 in overlaps[s1]:
                        ov = model.addVar(vtype=GRB.BINARY, name=f"Overlap3000_{c1}_{c2}_{s1}_{s2}")
                        model.addConstr(ov >= x[c1, s1] + x[c2, s2] - 1)

    for c_design in design_courses:
        if c_design in courses:
            for c_3000 in courses_3000:
                if c_3000 in courses:
                    for s1 in slots:
                        if s1 in overlaps:
                            for s2 in overlaps[s1]:
                                ov = model.addVar(vtype=GRB.BINARY, name=f"OverlapDesign_{c_design}_{c_3000}_{s1}_{s2}")
                                model.addConstr(ov >= x[c_design, s1] + x[c_3000, s2] - 1)
            for c_cap in df_capstone["course"].tolist():
                if c_cap in courses:
                    for s1 in slots:
                        if s1 in overlaps:
                            for s2 in overlaps[s1]:
                                ov = model.addVar(vtype=GRB.BINARY, name=f"OverlapDesignCapstone_{c_design}_{c_cap}_{s1}_{s2}")
                                model.addConstr(ov >= x[c_design, s1] + x[c_cap, s2] - 1)

    for c_tech in tech_electives:
        if c_tech in courses:
            for c_cap in df_capstone["course"].tolist():
                if c_cap in courses:
                    for s1 in slots:
                        if s1 in overlaps:
                            for s2 in overlaps[s1]:
                                ov = model.addVar(vtype=GRB.BINARY, name=f"OverlapTech_{c_tech}_{c_cap}_{s1}_{s2}")
                                model.addConstr(ov >= x[c_tech, s1] + x[c_cap, s2] - 1)

    # Objective
    penalty_slot = -1
    penalty_overlap_3000 = -1
    penalty_overlap_design = -1
    penaly_overlap_tech = -1

    model.setObjective(
        gp.quicksum(z[c] for c in courses) +
        gp.quicksum(u_var for u_var in u.values()) +
        penalty_slot * gp.quicksum(w_penalty[c, s] for c in courses for s in penalty_slots) +
        penalty_overlap_3000 * gp.quicksum([v for v in model.getVars() if v.VarName.startswith("Overlap3000_")]) +
        penalty_overlap_design * gp.quicksum([v for v in model.getVars() if v.VarName.startswith("OverlapDesign_")]) +
        penaly_overlap_tech * gp.quicksum([v for v in model.getVars() if v.VarName.startswith("OverlapTech_")]),
        GRB.MAXIMIZE
    )

    # Optimize
    model.optimize()

    # Read back Gurobi log (for non-streaming endpoint or debugging)
    log_text = ""
    try:
        with open(LOG_PATH, "r", encoding="utf-8", errors="ignore") as f:
            log_text = f.read()
        if len(log_text) > 200_000:
            log_text = log_text[-200_000:]
    except Exception:
        log_text = ""

    # Build solution
    if model.SolCount == 0:
        return pd.DataFrame(columns=["Course", "Day", "Start", "End"]), log_text

    model.Params.SolutionNumber = 0

    slots_list = df_slots["SlotID"].tolist()
    schedule = []
    for course in courses:
        assigned_slots = [s for s in slots_list if (model.getVarByName(f"x[{course},{s}]").Xn > 0.5)]
        for s in assigned_slots:
            row = df_slots.loc[df_slots["SlotID"] == s]
            schedule.append({
                "Course": course,
                "Day": row["Day"].values[0],
                "Start": row["Start"].values[0],
                "End": row["End"].values[0],
            })
    df_schedule = pd.DataFrame(schedule)

    if df_schedule.empty:
        return pd.DataFrame(columns=["Course", "Day", "Start", "End"]), log_text

    df_day_combined = (
        df_schedule.groupby(["Course", "Day"], as_index=False)
        .agg({"Start": "min", "End": "max"})
    )

    day_numbers = {"Monday": 1, "Tuesday": 2, "Wednesday": 3, "Thursday": 4, "Friday": 5}
    day_short = {"Monday": "Mo", "Tuesday": "Tu", "Wednesday": "We", "Thursday": "Th", "Friday": "Fr"}

    def _aggregate(group):
        return pd.Series({
            "Day": ",".join([day_short[d] for d in sorted(group["Day"], key=lambda d: day_numbers[d])]),
            "Start": group["Start"].min(),
            "End": group["End"].max()
        })

    df_final = (
        df_day_combined.groupby("Course", as_index=False)[["Day", "Start", "End"]]
        .apply(_aggregate)
        .reset_index(drop=True)
    )
    df_final = df_final[["Course", "Day", "Start", "End"]]
    return df_final, log_text


def run_and_optionally_write_excel(
    plan_semester: str,
    offerings_path: str,
    input_data_path: str,
    output_excel_path: Optional[str]
):
    df_final, log_text = run_optimization(plan_semester, offerings_path, input_data_path)

    # build Excel in-memory; write if path provided
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
        df_final.to_excel(writer, sheet_name="Optimal_Solution", index=False)
    excel_bytes = buffer.getvalue()

    if output_excel_path:
        with open(output_excel_path, "wb") as f:
            f.write(excel_bytes)

    return df_final, excel_bytes, log_text
