# -*- coding: utf-8 -*-
import pandas as pd
import re

# Load the Excel file of previous offerings
prev_offerings = pd.ExcelFile('Hx_Data_Offerings.xlsx')

# Extract the order of semesters dynamically from sheet names
semester_order = prev_offerings.sheet_names

# Create an empty list to gather data from all sheets
all_data = []

# Process each sheet to collect Meeting Pattern data
for sheet in semester_order:
    df = prev_offerings.parse(sheet)
    if df.empty:
        continue  # Skip empty sheets
    df['Semester'] = sheet  # Add a column to keep track of the semester
    # Combine 'Subject Code' and 'Catalog Number' to create 'Course'
    df['Course'] = df['Subject Code'] + " " + df['Catalog Number'].astype(str)
    all_data.append(df[['Course', 'Meeting Pattern', 'Semester']])

# Combine all sheets into a single DataFrame
combined_data = pd.concat(all_data, ignore_index=True)

# Sort by Semester in reverse order of recency
combined_data['Semester'] = pd.Categorical(combined_data['Semester'], categories=semester_order, ordered=True)
combined_data = combined_data.sort_values(by=['Course', 'Semester'], ascending=[True, False])

# Ensure that preference is assigned when both "Does Not Meet" and a valid meeting pattern exist
combined_data = combined_data.sort_values(by=['Course', 'Meeting Pattern'], key=lambda x: x != "Does Not Meet").drop_duplicates(subset=['Course'], keep='last')

# Extract the latest Meeting Pattern for each course
latest_data = combined_data.groupby('Course').first().reset_index()

# Define a function to extract the first valid duration
def extract_first_duration(meeting_pattern):
    if pd.isna(meeting_pattern) or "Does Not Meet" in meeting_pattern:
        return None
    try:
        # Use regex to extract the first valid time range from the string
        time_match = re.search(r'(\d{1,2}(:\d{2})?[ap]m)-(\d{1,2}(:\d{2})?[ap]m)', meeting_pattern)
        if time_match:
            start_time, _, end_time, _ = time_match.groups()
            # Add ":00" if minutes are missing
            if ':' not in start_time:
                start_time = start_time.replace("am", ":00am").replace("pm", ":00pm")
            if ':' not in end_time:
                end_time = end_time.replace("am", ":00am").replace("pm", ":00pm")
            # Convert to datetime for calculation
            start_time = pd.to_datetime(start_time.strip(), format='%I:%M%p')
            end_time = pd.to_datetime(end_time.strip(), format='%I:%M%p')
            # Calculate the duration in minutes
            duration = (end_time - start_time).seconds // 60
            return duration
    except Exception as e:
        print(f"Error processing pattern '{meeting_pattern}': {e}")
        return None

# Apply the function to calculate durations
latest_data['Duration'] = latest_data['Meeting Pattern'].apply(extract_first_duration)

# Determine preferences for 50-minute and 80-minute durations
latest_data['50_min'] = (latest_data['Duration'] == 50).astype(int)
latest_data['80_min'] = (latest_data['Duration'] == 80).astype(int)

# Prepare the final dataframe with Course, 50_min, and 80_min columns
df_preference = latest_data[['Course', '50_min', '80_min']]
df_preference.to_excel('preference.xlsx', index=False)

