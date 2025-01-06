using Pkg
Pkg.activate(".")
Pkg.instantiate() # Shift + Enter to execute a single line of code 
Pkg.status()
# Highlight the lines you want to execute and press Shift + Enter to execute the multiple lines

using CSV
using DataFrames
using Dates
using Statistics
using PlotlyJS
using Dash

# Reading CSV files into a DataFrame
df = CSV.read("data/SOFI_all_data.csv", DataFrame); # You can leave comments by typing the hash key
df2 = CSV.read("data/SOFI_coordinates.csv", DataFrame); # The semicolon ; at the end of a line suppresses the output of the command

# Sorting the DataFrame based on the values in the column named :timestamp
sort!(df, :timestamp)

available_indicators = sort(unique(df.station_id))

# Mapping station IDs to names
station_names = Dict(
    "SOFIS1_A8610A3135399218" => "SOFIS1",
    "SOFIS2_A8610A3135428119" => "SOFIS2",
    "SOFIS3_A8610A3135197407" => "SOFIS3",
    "SOFIS4_A8610A3439296D08" => "SOFIS4"
)


include("config.jl") 

# Maping of the coordinates
lon = df2.longitude_deg
lat = df2.latitude_deg
key = [station_names[k] for k in df2.station_id]

trace = scattermapbox(
    lat = lat,
    lon = lon,
    customdata = key,
    marker_size = 10,
    colorscale = "YlOrRd",
    hovertemplate = """
                      <b>%{customdata}</b><br>
                      Latitude: %{lat}<br>
                      Longitude: %{lon}
                      <extra></extra>
                      """
)

# Visualization settings
layout = Layout(
    height = "500",
    mapbox_accesstoken = MAPBOX_TOKEN, # Use your own Mapbox Access Token here (see README.md)
    mapbox_center_lon = median(lon),
    mapbox_center_lat = median(lat),
    mapbox_zoom = 13,
    mapbox_style = "mapbox://styles/mapbox/satellite-v9"  # Set Satellite as the default style
)

# Dropdown menu for map styles
dropdown_buttons = [
    Dict(
        "args" => ["mapbox.style", "mapbox://styles/mapbox/satellite-v9"],
        "label" => "Satellite",
        "method" => "relayout"
    ),
    Dict(
        "args" => ["mapbox.style", "mapbox://styles/mapbox/streets-v12"],
        "label" => "Streets",
        "method" => "relayout"
    )
]

updatemenus = [
    Dict(
        "buttons" => dropdown_buttons,
        "direction" => "down",
        "showactive" => true,
        "x" => 0.02,  # Adjusted x coordinate
        "xanchor" => "left",
        "y" => 0.95,  # Adjusted y coordinate
        "yanchor" => "top"
    )
]

layout["updatemenus"] = updatemenus

# Display figure
fig = plot(trace, layout)


# Convert timestamp to Date
df.date = Dates.Date.(df.timestamp)

# Filter data from April 15th to September 24th
df = filter(row -> 4 ≤ Dates.month(row.date) ≤ 9 &&
                    !((Dates.month(row.date) == 4 && Dates.day(row.date) < 15) ||
                    (Dates.month(row.date) == 9 && Dates.day(row.date) > 24)), df)

# Group by date and station_id and aggregate the data
daily_data = combine(groupby(df, [:date, :station_id])) do grouped
    DataFrame(
        date = first(grouped.date),
        station_id = first(grouped.station_id),
        precipitation_mm = round(sum(grouped.precipitation_accum_mm), digits=2),
        avg_temperature_c = round(mean(grouped.temperature_c), digits=2),
        min_temperature_c = round(minimum(grouped.temperature_c), digits=2),
        max_temperature_c = round(maximum(grouped.temperature_c), digits=2),
        avg_humidity_pct = round(mean(grouped.rel_humidity_pctg), digits=2),
        avg_wind_speed_kmh = round(mean(grouped.wind_speed_kmh), digits=2)
    )
end

# 1 scenario 10 mm precipitation 
daily_data[!, :occurrence_1] .= 0
for i in 1:size(daily_data, 1) - 1
    if daily_data[i, :precipitation_mm] + daily_data[i + 1, :precipitation_mm] ≥ 10 &&
        daily_data[i + 1, :min_temperature_c] ≥ 10
        daily_data[i + 1, :occurrence_1] = 1
    else
        daily_data[i + 1, :occurrence_1] = 0
    end
end

# 2 scenario 7 mm precipitation 
daily_data[!, :occurrence_2] .= 0
for i in 1:size(daily_data, 1) - 1
    if daily_data[i, :precipitation_mm] + daily_data[i + 1, :precipitation_mm] ≥ 7 &&
        daily_data[i + 1, :min_temperature_c] ≥ 10
        daily_data[i + 1, :occurrence_2] = 1
    else
        daily_data[i + 1, :occurrence_2] = 0
    end
end

# 3 scenario 5 mm precipitation
daily_data[!, :occurrence_3] .= 0
for i in 1:size(daily_data, 1) - 1
    if daily_data[i, :precipitation_mm] + daily_data[i + 1, :precipitation_mm] ≥ 5 &&
        daily_data[i + 1, :min_temperature_c] ≥ 10
        daily_data[i + 1, :occurrence_3] = 1
    else
        daily_data[i + 1, :occurrence_3] = 0
    end
end

# External stylesheet for Bootstrap
external_stylesheets = ["https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css"]

# Initialize the Dash app with external stylesheets
app = dash(external_stylesheets=external_stylesheets)

# Define the layout of the app
app.layout = html_div() do
    html_div(style=Dict("display" => "flex", "align-items" => "center", "justify-content" => "flex-end", "margin" => "20px"), children=[
    ]),
    html_div(style=Dict("text-align" => "center", "margin" => "20px"), children=[
        html_h1("Prototype Vineyard Monitoring Dashboard"),
        html_h3("Vineyard Monitoring Made Simple 🍇", style=Dict("color" => "grey")),
    ]),
    html_h4("Location of the installed SOFIS sensors", style=Dict("margin" => "40px")),
    dcc_graph(id="density-map", figure=fig),
    html_h4("Weather data", style=Dict("margin" => "40px")),
    html_div(children=[
        html_div(
            children=[
                dcc_dropdown(
                    id="sensor-dropdown",
                    options=[Dict("label" => station_names[i], "value" => i) for i in available_indicators],
                    value="SOFIS1_A8610A3135399218",
                    style=Dict("width" => "100%")
                ),
            ],
            style=Dict("width" => "48%", "display" => "inline-block", "margin" => "40px")
        ),
        dcc_graph(id="temperature-graphic", style=Dict("margin" => "20px")),
        dcc_graph(id="humidity-graphic", style=Dict("margin" => "20px")),
        dcc_graph(id="precipitation-graphic", style=Dict("margin" => "20px")),
        dcc_graph(id="radiation-graphic", style=Dict("margin" => "20px")),
        dcc_graph(id="wind-speed-graphic", style=Dict("margin" => "20px")),
    ]),
    html_h4("Downy mildew", style=Dict("margin" => "40px")),
    dcc_markdown("""
        Downy mildew, caused by *Plasmopara viticola*, is a major disease of grapevines in the Azores. It attacks all green parts of the grapevine, including branches, leaves, clusters, and tendrils. Typically, its presence is first noticed on the leaves, characterized by the appearance of oil spots. In severe infestations, these spots may progress to partial or complete desiccation. Additionally, symptoms can manifest in the rachis, which darkens and curves into an 'S' shape, eventually becoming coated with whitish dust. Similar symptoms can also occur on inflorescences and newly formed grapes.

        Downy mildew is driven by the weather. The 10s method helps predict the primary mildew infection of grapevines.

        The **10s method** refers to **10 cm vine shoot**, **10ºC temperature**, and **10 mm of accumulated rainfall** within two consecutive days.

        Due to the unique microclimate in the Azores and the uncertainty regarding fungal diseases, we also test *7 mm* and *5 mm* rainfall scenarios. This allows producers to choose the scenario that best suits their conditions.
    """, style=Dict("margin" => "40px", "text-align" => "justify")),

    html_div(
        children=[
            html_div(
                children=[
                    dcc_dropdown(
                        id="sensor-dropdown2",
                        options=[Dict("label" => station_names[station], "value" => station) for station in unique(daily_data.station_id)],
                        value=unique(daily_data.station_id)[1],  # Default value
                        style=Dict("width" => "100%")
                    ),
                ],
                style=Dict("width" => "48%", "display" => "inline-block", "margin" => "40px")
            ),
            dcc_graph(id="primary-infection-graph", style=Dict("margin" => "20px")),
        ],
        style=Dict("width" => "100%", "display" => "inline-block")
    ),
    dcc_markdown("""
        **Note that the model should only be considered if the vine shoot is 10 cm and data is available.** 
        Measure the vine shoots regularly and ensure the sensor is turned on and functioning properly, i.e., 
        there is at least temperature data available in the weather data section. The model applies 
        from April 15th to September 15th, when the vine has green parts.
    """, style=Dict("margin" => "40px", "text-align" => "justify")),
    html_div(
        style=Dict("display" => "flex", "align-items" => "center", "justify-content" => "center"),
        children=[
            html_img(src="assets/aircentre_logoTAG_rgb-01.png", style=Dict("width" => "100px", "height" => "auto")),
            html_img(src="assets/logo-cmah.png", style=Dict("width" => "100px", "height" => "auto", "margin-left" => "30px", "margin-right" => "30px")),
            html_img(src="assets/sfcolab.png", style=Dict("width" => "100px", "height" => "auto"))
        ]
    ),
    html_p("Copyright © 2024 AIR Centre. All rights reserved.", style=Dict("text-align" => "center", "margin" => "20px"))
end

# Callbacks
callback!(
    app,
    Output("temperature-graphic", "figure"),
    Input("sensor-dropdown", "value")
) do sensor_value
    dff_sensor = df[df.station_id .== sensor_value, :]
    temp_plot = Plot(
        scatter(
            x=dff_sensor.timestamp,
            y=dff_sensor.temperature_c,
            mode="markers",
            marker=attr(size=5, opacity=0.5, line=attr(width=0.5, color="white"))
        ),
        Layout(
            title=attr(
                text="Temperature",
                x=0.5,
                xanchor="center",
                yanchor="top",
                pad=attr(t=30)  # Add padding above the title
            ),
            xaxis=attr(
                rangeselector=attr(
                    buttons=[
                        attr(count=1, label="1M", step="month", stepmode="backward"),
                        attr(count=3, label="3M", step="month", stepmode="backward"),
                        attr(count=6, label="6M", step="month", stepmode="backward"),
                        attr(count=1, label="YTD", step="year", stepmode="todate"),
                        attr(count=1, label="1Y", step="year", stepmode="backward"),
                        attr(label="All", step="all")
                    ],
                    pad=attr(t=30)  # Add padding above the range selector
                )
            ),
            yaxis=attr(
                title="Temperature (°C)"
            ),
            hovermode="closest",
            margin=attr(t=100)  # Adjust top margin to ensure space for the title and buttons
        )
    )
    return temp_plot
end

callback!(
    app,
    Output("humidity-graphic", "figure"),
    Input("sensor-dropdown", "value")
) do sensor_value
    dff_sensor = df[df.station_id .== sensor_value, :]
    humidity_plot = Plot(
        scatter(
            x=dff_sensor.timestamp,
            y=dff_sensor.rel_humidity_pctg,
            mode="markers",
            marker=attr(size=5, opacity=0.5, line=attr(width=0.5, color="white"))
        ),
        Layout(
            title=attr(
                text="Relative Humidity",
                x=0.5,
                xanchor="center",
                yanchor="top",
                pad=attr(t=30)
            ),
            xaxis=attr(
                rangeselector=attr(
                    buttons=[
                        attr(count=1, label="1M", step="month", stepmode="backward"),
                        attr(count=3, label="3M", step="month", stepmode="backward"),
                        attr(count=6, label="6M", step="month", stepmode="backward"),
                        attr(count=1, label="YTD", step="year", stepmode="todate"),
                        attr(count=1, label="1Y", step="year", stepmode="backward"),
                        attr(label="All", step="all")
                    ],
                    pad=attr(t=30)
                )
            ),
            yaxis=attr(
                title="Relative Humidity (%)"
            ),
            hovermode="closest",
            margin=attr(t=100)
        )
    )
    return humidity_plot
end

callback!(
    app,
    Output("precipitation-graphic", "figure"),
    Input("sensor-dropdown", "value")
) do sensor_value
    dff_sensor = df[df.station_id .== sensor_value, :]
    precipitation_plot = Plot(
        bar(
            x=dff_sensor.timestamp,
            y=dff_sensor.precipitation_accum_mm,
            marker=attr(color="blue", opacity=1.0)
        ),
        Layout(
            title=attr(
                text="Precipitation",
                x=0.5,
                xanchor="center",
                yanchor="top",
                pad=attr(t=30)
            ),
            xaxis=attr(
                rangeselector=attr(
                    buttons=[
                        attr(count=1, label="1M", step="month", stepmode="backward"),
                        attr(count=3, label="3M", step="month", stepmode="backward"),
                        attr(count=6, label="6M", step="month", stepmode="backward"),
                        attr(count=1, label="YTD", step="year", stepmode="todate"),
                        attr(count=1, label="1Y", step="year", stepmode="backward"),
                        attr(label="All", step="all")
                    ],
                    pad=attr(t=30)
                )
            ),
            yaxis=attr(
                title="Precipitation (mm)"
            ),
            hovermode="closest",
            margin=attr(t=100)
        )
    )
    return precipitation_plot
end

callback!(
    app,
    Output("radiation-graphic", "figure"),
    Input("sensor-dropdown", "value")
) do sensor_value
    dff_sensor = df[df.station_id .== sensor_value, :]
    radiation_plot = Plot(
        bar(
            x=dff_sensor.timestamp,
            y=dff_sensor.radiation_kjm2,
            marker=attr(color="blue", opacity=1.0)
        ),
        Layout(
            title=attr(
                text="Photosynthetically Active Radiation (PAR)",
                x=0.5,
                xanchor="center",
                yanchor="top",
                pad=attr(t=30)
            ),
            xaxis=attr(
                rangeselector=attr(
                    buttons=[
                        attr(count=1, label="1M", step="month", stepmode="backward"),
                        attr(count=3, label="3M", step="month", stepmode="backward"),
                        attr(count=6, label="6M", step="month", stepmode="backward"),
                        attr(count=1, label="YTD", step="year", stepmode="todate"),
                        attr(count=1, label="1Y", step="year", stepmode="backward"),
                        attr(label="All", step="all")
                    ],
                    pad=attr(t=30)
                )
            ),
            yaxis=attr(
                title="PAR (mmol quantum / m²⋅s)"
            ),
            hovermode="closest",
            margin=attr(t=100)
        )
    )
    return radiation_plot
end

callback!(
    app,
    Output("wind-speed-graphic", "figure"),
    Input("sensor-dropdown", "value")
) do sensor_value
    dff_sensor = df[df.station_id .== sensor_value, :]
    wind_speed_plot = Plot(
        bar(
            x=dff_sensor.timestamp,
            y=dff_sensor.wind_speed_kmh,
            marker=attr(color="blue", opacity=1.0)
        ),
        Layout(
            title=attr(
                text="Wind Speed",
                x=0.5,
                xanchor="center",
                yanchor="top",
                pad=attr(t=30)
            ),
            xaxis=attr(
                rangeselector=attr(
                    buttons=[
                        attr(count=1, label="1M", step="month", stepmode="backward"),
                        attr(count=3, label="3M", step="month", stepmode="backward"),
                        attr(count=6, label="6M", step="month", stepmode="backward"),
                        attr(count=1, label="YTD", step="year", stepmode="todate"),
                        attr(count=1, label="1Y", step="year", stepmode="backward"),
                        attr(label="All", step="all")
                    ],
                    pad=attr(t=30)
                )
            ),
            yaxis=attr(
                title="Wind Speed (km/h)"
            ),
            hovermode="closest",
            margin=attr(t=100)
        )
    )
    return wind_speed_plot
end

callback!(
    app,
    Output("primary-infection-graph", "figure"),
    Input("sensor-dropdown2", "value")
) do selected_sensor
    sensor_data = filter(row -> row.station_id == selected_sensor, daily_data)
    sort!(sensor_data, :date)

    # Extract necessary data
    dates = sensor_data.date
    min_temperature = sensor_data.min_temperature_c
    precipitation = sensor_data.precipitation_mm

    # Create traces
    trace_temp = scatter(
        x=dates,
        y=min_temperature,
        mode="lines+markers",
        name="Daily Minimum Temperature (°C)",
        line=attr(color="green")
    )
    trace_precipitation = bar(
        x=dates,
        y=precipitation,
        name="Daily Accumulated Precipitation (mm)",
        marker=attr(color="blue")
    )
    threshold_line = scatter(
        x=dates,
        y=fill(10, length(dates)),
        mode="lines",
        line=attr(color="red"),
        name="10°C Threshold"
    )

    # Occurrence markers
    occurrence_1_dates = sensor_data[sensor_data.occurrence_1 .== 1, :date]
    trace_occurrence_1 = scatter(
        x=occurrence_1_dates,
        y=fill(10, length(occurrence_1_dates)),
        mode="markers",
        name="Threshold (10 mm)",
        marker=attr(color="orange")
    )

    occurrence_2_dates = sensor_data[sensor_data.occurrence_2 .== 1, :date]
    trace_occurrence_2 = scatter(
        x=occurrence_2_dates,
        y=fill(7, length(occurrence_2_dates)),
        mode="markers",
        name="Threshold (7 mm)",
        marker=attr(color="purple")
    )

    occurrence_3_dates = sensor_data[sensor_data.occurrence_3 .== 1, :date]
    trace_occurrence_3 = scatter(
        x=occurrence_3_dates,
        y=fill(5, length(occurrence_3_dates)),
        mode="markers",
        name="Threshold (5 mm)",
        marker=attr(color="yellow")
    )

    # Combine all traces
    plot_data_combined = [
        trace_temp,
        trace_precipitation,
        threshold_line,
        trace_occurrence_1,
        trace_occurrence_2,
        trace_occurrence_3
    ]

    # Layout
    layout_combined = Layout(
        title=attr(
            text="Primary Infection Model of Downy Mildew",
            x=0.5,
            xanchor="center",
            yanchor="top",
            pad=attr(t=30)
        ),
        yaxis=attr(
            title="Value"
        ),
        hovermode="closest",
        margin=attr(t=100)
    )
    
    # Create the figure
    fig_combined = Plot(plot_data_combined, layout_combined)
    return fig_combined
end

run_server(app, "0.0.0.0", 8050, debug=true)

# Open your web browser and go to http://127.0.0.1:8050/ to view the dashboard locally

# Ctrl + C to stop the dash app

