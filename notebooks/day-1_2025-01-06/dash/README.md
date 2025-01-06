# Real-time alert systems in Julia

**Author:** [Iga Szczesniak](https://igaszczesniak.github.io)  
**Location:** Terceira Island, Azores, Portugal <br>
**Date:** January 6th 2025


This repository provides resources for creating interactive, web-based data applications in Julia using  Dash.jl framework. The example project featured here is the [**vineyard monitoring dashboard**](https://services.aircentre.org/agrodigital/terceira/) (but in English), which includes interactive elements such as dropdown lists, maps, graphs, and buttons.

## Prerequisites

Ensure you have the following software installed before setting up the project:

- [Visual Studio Code (VS Code)](https://code.visualstudio.com/download)
- [Julia v1.11](https://julialang.org/downloads/#upcoming_release)

## Setup Instructions

### 1. Clone the Repository

Open your terminal and clone the repository by running:

```bash
git clone https://github.com/igaszczesniak/dashboard-julia.git
```
### 2. Mapbox Access Token Configuration

Mapbox is used for styling maps within the dashboard. To configure it, follow these steps:

#### Acquire a Mapbox Access Token

1. Visit the [Mapbox Access Tokens page](https://account.mapbox.com/access-tokens/clv3yseh502zv2in19d765bww/).
2. Log in or sign up for a Mapbox account.
3. Generate a new access token or use an existing one.

#### Configure the Token in Your Project

1. Navigate to the project folder.
2. Create a file named `config.jl` in the project folder.
3. Open `config.jl` and define your Mapbox Access Token by adding the following line:

    ```julia
    MAPBOX_TOKEN = "paste-your-access-token-here"
    ```

    Replace `"paste-your-access-token-here"` with the actual token you obtained from Mapbox.

## Running the Project

### 1. Launch the Dash App

Open the project in VS Code and run the following command in your terminal to start the Dash app:

```julia
julia dashboard.jl
```

### 2. View the Dashboard

Open your web browser and navigate to http://127.0.0.1:8050/ to view the dashboard locally!