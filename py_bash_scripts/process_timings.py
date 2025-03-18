import glob
import re
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

data = []

# Regex pattern: one or more digits optionally followed by a dot and more digits.
num_pattern = r'\d+(?:\.\d+)?'

# Process each file matching the pattern.
for filepath in glob.glob("timings_*.txt"):
    # Extract the file tag: the part after the first underscore and before the dot.
    file_tag_match = re.search(r'timings_(.*?)\.txt', filepath)
    file_tag = file_tag_match.group(1) if file_tag_match else "unknown"
    
    with open(filepath, "r") as f:
        lines = f.readlines()
    
    # Iterate over lines looking for consecutive pair of lines.
    i = 0
    while i < len(lines):
        if "Total new embeddings processed in this run:" in lines[i]:
            # Extract the embedding number using the stricter regex.
            embeddings_nums = re.findall(num_pattern, lines[i])
            if i + 1 < len(lines) and "Total time:" in lines[i+1]:
                # Extract numbers from the "Total time:" line.
                time_nums = re.findall(num_pattern, lines[i+1])
                # In case there are more than three numbers, select only the first three.
                time_nums = time_nums[:3]
                # Combine the file tag, embeddings number, and time numbers.
                data.append([file_tag] + embeddings_nums + time_nums)
                i += 2  # Skip the next line as it is already processed.
                continue
        i += 1

# Create a DataFrame with one row per processed pair.
df = pd.DataFrame(data, columns=["file_tag", "embeddings", "total_time", "time_per_prot", "avg_len"])

print(df)

convert the embeddings column to numeric and filter out rows where it equals 0.
df["embeddings"] = pd.to_numeric(df["embeddings"])
df_filtered = df[df["embeddings"] != 0]

# Write the filtered DataFrame to CSV.
csv_filename = "processed_timings.csv"
df_filtered.to_csv(csv_filename, index=False)
print(f"CSV file written to: {csv_filename}")

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# === 1. Load Data and Order file_tag ===
# Assume df_read is the DataFrame read from processed_timings.csv.
df_read = pd.read_csv("processed_timings.csv")

# Define the natural ordering of file_tags.
cat_order = [
    "0_99", "100_199", "200_299", "300_399", "400_499",
    "500_599", "600_699", "700_799", "800_899", "900_1000", "1001min"
]

# Set file_tag as a categorical with the desired order.
df_read["file_tag"] = pd.Categorical(df_read["file_tag"], categories=cat_order, ordered=True)

# === 2. Compute Summary Statistics for time_per_prot grouped by file_tag ===
summary_stats = (
    df_read
    .groupby("file_tag", observed=True)["time_per_prot"]
    .agg(["min", "median", "max", "mean", "std"])
    .reset_index()
)

# === 3. Fit a Quadratic Curve for the Scatter Plot ===
# Use avg_len as x and time_per_prot as y.
x_data = df_read["avg_len"].values
y_data = df_read["time_per_prot"].values

# Fit a 2nd-degree polynomial: y = a*x² + b*x + c
coeffs = np.polyfit(x_data, y_data, deg=2)  # [a, b, c]
poly = np.poly1d(coeffs)
a, b, c = coeffs

# Generate fitted curve data
x_lin = np.linspace(x_data.min(), x_data.max(), 2*(int(x_data.max()) - int(x_data.min())))
y_lin = poly(x_lin)

# Create multiline text for the equation (using HTML <br> for line breaks)
equation_text = f"y = {a:.6g}x²<br>+ {b:.6g}x<br>+ {c:.6g}"

# === 4. Create the Scatter Plot with Plotly Express and Add the Fit Line ===
fig = px.scatter(
    df_read,
    x="avg_len",
    y="time_per_prot",
    color="file_tag",
    category_orders={"file_tag": cat_order},  # Force legend order
    hover_data=["embeddings", "total_time"],
    title="Time per Protein vs. Average Protein Length with Fit and Summary Stats",
    labels={
        "avg_len": "Average Protein Length",
        "time_per_prot": "Time per Protein (s)",
        "file_tag": "File Tag"
    }
)

# Add the quadratic fit line
fig.add_trace(
    go.Scatter(
        x=x_lin,
        y=y_lin,
        mode='lines',
        name='Quadratic Fit',
        line=dict(color='black', dash='dash')
    )
)

# Get the quadratic fit trace (assuming it's the last one added)
qfit_trace = fig.data[-1]
# Get all the other scatter traces
scatter_traces = list(fig.data[:-1])
# Reassemble fig.data with the quadratic fit trace first
fig.data = (qfit_trace,) + tuple(scatter_traces)

# Add the equation annotation in the bottom-right with a bounding box
fig.add_annotation(
    xref="paper", yref="paper",
    x=0.9, y=0.1,  # Position near bottom-right
    text=equation_text,
    showarrow=False,
    font=dict(size=12),
    bordercolor="black",
    borderwidth=1,
    borderpad=5,
    bgcolor="white",
    align="left"
)

# === 5. Add the Summary Stats as an Inset Table in the Top-Left ===
formatted_values = []
for col in summary_stats.columns:
    if pd.api.types.is_numeric_dtype(summary_stats[col]):
        formatted_values.append(summary_stats[col].apply(lambda x: f"{x:.3g}").tolist())
    else:
        formatted_values.append(summary_stats[col].tolist())

# Create a table trace from the summary_stats DataFrame.
table_trace = go.Table(
    header=dict(
        values=list(summary_stats.columns),
        fill_color='paleturquoise',
        align='left'
    ),
    cells=dict(
        values=formatted_values,
        fill_color='lavender',
        align='left'
    ),
    # Position the table in an inset fashion (domain in paper coordinates)
    domain=dict(x=[0.01, 0.4], y=[0.4, 0.97])
)

# Add the table trace to the figure.
fig.add_trace(table_trace)

# === 6. Export the Plot as an Interactive HTML File ===
html_filename = "length_vs_time_interactive_plot.html"
fig.write_html(html_filename)
print(f"Plot exported to: {html_filename}")

# Show the final plot
fig.show()

