### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 97a4183e-f670-47b6-a318-4254f34cac32
begin
	using JLSO
	using Flux
	using Distributions
	using PlutoPlotly
	using BenchmarkTools
	using DataFrames
end

# ╔═╡ a9c6cecd-0013-42b5-b7cc-0362bdb4a663
md"""
# Introduction
"""

# ╔═╡ 7e34c342-91dd-4c15-8a4c-0d6a61093371
md"""
# Loading Packages
- `JLSO` will be used to deserialise data to train the surrogate
- `Flux` is a ML library that will let us train the network and do the downstream task
- `Distributions` give us functionality to create nice initializations for our model
- `PlutoPlotly` for beautiful plots to visualize and debug our surrogate's performance
- `BenchmarkTools` to benchmark the speed of the surrogate and simulation
- `DataFrames` used to conveniently organise data for the parallel coordinate plot
"""

# ╔═╡ 62dba600-c90b-11ef-009f-07f05e9e18cc
md"""
#  Loading Data

Data was generated using RainMaker.jl and SpeedyWeather.jl. The procedure can be summarised as follows:

1. Specify upper and lower bounds, including number of samples
2. Use a Quasi Monte Carlo scheme such as `LatinHypercubeSample` to sample this space
3. Simulate each sample in parallel and return a tuple that is (params,output)
4. Serialise the data for safe keeping (optional)

The code used has been provided as markdown for convenience but for this workshop, we are going to be using a dataset already generated in the interest of time.

```julia
using Distributed
addprocs(64)
using QuasiMonteCarlo
@everywhere using RainMaker    
@everywhere using SpeedyWeather
include("utils.jl") #where max_precipitation lives


@everywhere const PARAMETER_KEYS = (
    :orography_scale,           # [1],      default: 1, scale of global orography
    :mountain_height,           # [m],      default: 0, height of an additional azores mountain
    :mountain_size,             # [˚],      default: 1, horizontal size of an additional azores mountain
    :mountain_lon,              # [˚E],     default: -27.25, longitude of an additional azores mountain
    :mountain_lat,              # [˚N],     default: 38.7, latitude of an additional azores mountain
    :temperature_equator,       # [K],      default: 300, sea surface temperature at the equator
    :temperature_pole,          # [K],      default: 273, sea surfaec temperature at the poles
    :temperature_atlantic,      # [K],      default: 0, sea surface temperature anomaly in the atlantic
    :temperature_azores,        # [K],      default: 0, sea surface temperature anomaly at the azores
    :zonal_wind,                # [m/s],    default: 35, zonal wind speed
)

n = 10000
lb = [0, -2000, 0, -180, -90, 270, 270, -5, -5, 5]
ub = [2, 5000, 30, 180, 90, 300, 300, 5, 5, 50]

s = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample())

sols = pmap(eachcol(s)) do sample
    # Attempt the simulation
    sol = try
        max_precipitation(sample)
    catch e
        # If there's an error, store a fallback (e.g., NaN), 
        # or do some logging if you wish:
        @warn "Error in max_precipitation($sample): $e"
        NaN
    end
    return (sample, sol)
end

sampled_params = reduce(hcat,first.(sols))
sampled_outputs = last.(sols)

#save the data
using JLSO
JLSO.save("10kdata.jlso", Dict(:d=>(inputs = sampled_params, outputs = sampled_outputs)))
```
"""

# ╔═╡ a26dfb20-2644-4657-b5b5-864a32ac3053
begin
	data = JLSO.load("10kdata.jlso")[:d]
	input_data = data.inputs
	output_data = reshape(data.outputs, (1,:))
	nothing
end

# ╔═╡ af87db2b-31af-4aa9-b1a9-9d335bcde0e8
md"""
## Let's explore the 10k dataset
"""

# ╔═╡ 0c685bdf-edaf-4116-8f05-7e11383e3200
begin
	hist_plot = histogram(
	    x = vec(output_data),
	    name = "Distribution"
	)
	
	hist_layout = Layout(
	    title = attr(
	        text = "Total Precipitation - 10k LHC",
	        x = 0.5,             # Moves the title to the horizontal center
	        xanchor = "center"   
	    ),
	    xaxis = attr(title = "(mm)")
	)
	plot(hist_plot, hist_layout)
end

# ╔═╡ 08523ab9-5c66-4ceb-9c12-96ad1e5e9c01
extrema(output_data)

# ╔═╡ c13793fa-30ea-43e7-a416-7651877c140b
mean(output_data)

# ╔═╡ f1ee0a37-cae6-4223-a5b0-aea809fca768
var(output_data)

# ╔═╡ e57ee8db-fa0f-422f-8c0a-7ee612cc732a
md"""
## Discussion points:
- Data appears to have a lot of variability according to the mean and variance
- Histogram shows we have a heavy tail indicating we have behaviours that are under represented 
- Strong indications already that this dataset won't be the final one before we have a full fledged surrogate.
"""

# ╔═╡ 106e4cd3-e376-4c5b-b0eb-90952d79c760
md"""
# Preparing Data for Surrogate Training
1.  Normalize inputs and outputs to be between 0 and 1
2.  Split the data into `train` and `validation`
3.  Package data to be in batches and place on the gpu

"""

# ╔═╡ b5b681a6-f020-41c3-a690-51a9f730515b
md"""
# Initialise `surrogate` Model 
Choose Neural Network architecture (we defer to a simple Resnet for now) and specify hyperparameters (activation functions, hidden sizes etc..)

![Simple Resnet] (https://upload.wikimedia.org/wikipedia/commons/b/ba/ResBlock.png)

"""

# ╔═╡ 31e025b7-8b08-4adb-bc55-c4ebd72a931b
md"""
# Configure Training Hyperparameters
1. Use a basic schedule to decay learning rate over time
2. Decide on epochs per learning rate decrement
3. Choose the optimiser
4. track parameters for L2 regularization
"""

# ╔═╡ 22dad96e-c4c8-4ff8-a4d7-f83accb02b72
md"""
# Train the model!
1.  Run the model through the loop
2.  Sprinkle some standard logging (we look at the normalised training loss and the denormalised - original scale - validation loss)
"""

# ╔═╡ 70bfa951-27b3-4688-83c2-dd5be03ff244
xvals  = collect(0:100:100*500)[1:end-1] 

# ╔═╡ 607321b9-fd3f-41bf-91d0-2e7d099afc6c
md"""
### Parallel Coordinate Plot for Training Data
"""

# ╔═╡ df9b4742-f68e-4ff3-b501-ecbe3199d300
md"""
### Parallel Coordinate Plot for Valdation Data
"""

# ╔═╡ ca4456ab-dc1c-402c-98e1-07b7665bf86d
md"""
## Discussion Points
This particular part of the surrogate training shows how you can diagnose how well the model performs. Each line goes through a number for each parameter representing a sample. The end of the line is the mae for its prediction, and you can interact with the plot by dragging over the axis to select samples. Using this plot, one  can then perform active learning to generate more samples by examining where the model is performing poorly. We didn't train the model for that long so it's performance isn't very good. However, another model was trained in advanced that used a hidden size of 128, and trained for much longer which is provided in case compute isn't available for those trying this notebook and can't train something like that. As can be seen in it's parallel coordinate plots, we do much better on training and validation.
"""

# ╔═╡ 10e14cb4-ab21-4eda-b106-bc785c4ea859
pretrained_model = JLSO.load("128h-10k-2L-resnet.jlso")[:d]

# ╔═╡ 8b7d5e71-d10b-493e-addf-6df24971e385
md"""
### Parallel Coordinate Plot for Training Data: Pretrained
"""

# ╔═╡ 60770e3b-83a4-4a8a-a631-b92b0e0242d7
md"""
### Parallel Coordinate Plot for Validation Data: Pretrained
"""

# ╔═╡ 7c27f6c2-4869-4cee-9f49-b8530dde3c6f
md"""
# Downstream Optimization

Key idea here is to demonstrate the power of having a surrogate that is incredibly fast, is fully differentiable, and is batchable . Therefore we can do the following:

1. Randomly sample a series of parameters from the parameter space as starting points
2. Specify the objective of finding the parameters that result in the most total precipitation, and that the parameters stay in the bounds the network was trained on
3. Perform gradient descent to update the input parameters into the network that minimize this objective.
4. Choose from the sampled parameters that gave the largest output
5. Simulate it with the simulator to verify!
"""

# ╔═╡ f0441313-65d5-4a75-8499-2f17db7668a7
function input_objective(surr, X; λ=100.0)
	y_pred = surr(X) 
	obj = -sum(y_pred) #since we have to minimize the objective, the -sum will work here
	# Heavy differentiable penalty if X < 0 or X > 1
	lower_violation = @. max(0.0f0, -X)      
	upper_violation = @. max(0.0f0, X - 1.0f0) 
	penalty = sum(lower_violation.^2) + sum(upper_violation.^2)
	return obj + λ * penalty
end

# ╔═╡ 796e1216-fc68-41f0-b47b-1d638dff6bf4
begin
	X_candidate = rand(10,5)
	total_steps = 10000
	opt_in = Flux.Adam(1e-4)
	local st_opt = Flux.setup(opt_in, X_candidate)
	for step in 1:total_steps
	    gs = gradient(X_candidate) do guesses
	        input_objective(pretrained_model,guesses)
	    end
	    st_opt, X_candidate = Optimisers.update(st_opt, X_candidate, gs...)
	    # Print progress occasionally
	    if step % 200 == 0
	        @info "Step $step, objective() = $(input_objective(pretrained_model,X_candidate))"
	    end
	end
end

# ╔═╡ 199c8bf6-4ef7-43ee-8a44-d168dbf94304
md"""

### Verify!


```julia 
julia> max_precipitation(best_param_sample)
[ Info: RainGauge{Float32, AnvilInterpolator{Float32, OctahedralGaussianGrid}} callback added with key callback_xVEB
Weather is speedy: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████| Time: 0:00:04 (961.28 years/day)
229.42801f0
```

"""

# ╔═╡ 1f154133-62d4-47f6-8232-b5a35f927378
md"""
# Benchmarking Single Simulation vs `surrogate`

```julia
julia> @benchmark max_precipitation(best_inputs_orig)

BenchmarkTools.Trial: 2 samples with 1 evaluation.
 Range (min … max):  4.982 s …  4.990 s  ┊ GC (min … max): 0.93% … 1.52%
 Time  (median):     4.986 s             ┊ GC (median):    1.23%
 Time  (mean ± σ):   4.986 s ± 5.726 ms  ┊ GC (mean ± σ):  1.23% ± 0.42%

  █                                                      █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  4.98 s        Histogram: frequency by time        4.99 s <

 Memory estimate: 806.23 MiB, allocs estimate: 8323135.

julia> @benchmark surrogate(best_inputs_norm[:,1:1])

BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  19.671 μs … 447.771 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     24.760 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   25.069 μs ±   8.258 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

                              ▄▇█▄▂                              
  ▁▁▁▁▁▁▁▁▁▁▁▁▁▂▄▅▄▃▂▂▁▁▁▁▁▂▄▇█████▇▇▆▄▃▃▃▃▃▃▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  19.7 μs         Histogram: frequency by time         29.6 μs <

 Memory estimate: 7.69 KiB, allocs estimate: 12.
```
  
  With speedups like this you could probably sample 1000s or 10s of 1000s of parameters during the downstream optimization unlocking some serious throughput for the procedure 
  
"""



# ╔═╡ cae4b323-475f-44f4-8586-d3b41e1ae8b1
md"""
## Discussion points
With speedups like this you could probably sample 1000s or 10s of 1000s of parameters during the downstream optimization with half decent hardware unlocking some serious throughput for the procedure
"""

# ╔═╡ 29dd72ce-1145-436d-892c-c3b5d8693d69
md"""
Congratulations, you've gone through your first pass of generating and validating a surrogate! This is only the beginning of the journey, so many things to try to bring the surrogate's performance on the validation up to par with it's performance on training.

Things to try:
1. Sample much more data. This is the easiest and usually the first thing to try (think 100k - 1million samples)
2. Change the architecture. Add more layers, use different activation functions, play with regularisation
3. Sample data adaptively using active learning 


Feel free to try non-deep learning methods as well
- XGBoost.jl
- MLJ.jl
- LIBSVM.jl

"""

# ╔═╡ 5b4c9729-aafe-4231-9429-a11ada56baf1
md"""
# Convenience Functions
"""

# ╔═╡ 18aeca69-12d4-4414-932e-f87ba42f79a6
normalise(x,min,max) = @. (x-min) / (max-min)

# ╔═╡ 77e3c694-1731-4e14-98e8-b09cfbaa9f94
begin
# Normalise inputs (parameters), take bounds from earlier
	inputs_lb = [0, -2000, 0, -180, -90, 270, 270, -5, -5, 5]
	inputs_ub = [2, 5000, 30, 180, 90, 300, 300, 5, 5, 50]
	input_data_norm = normalise(input_data, inputs_lb, inputs_ub)
	nothing
end

# ╔═╡ 2c0b78e5-817f-4eea-b0e1-41d5cfcc8cc3
extrema(input_data_norm,dims=2)  #show input data is in between 0 and 1

# ╔═╡ f52bba40-fd7d-4646-82d3-cc6420219916
begin
	# Normalise outputs
	outputs_lb, outputs_ub = extrema(output_data)
	output_data_norm = normalise(output_data, outputs_lb, outputs_ub)
	nothing
end

# ╔═╡ 73ae3f0f-7bbf-44d5-8522-acce0c9afc3a
extrema(output_data_norm,dims=2)  #show output data is in between 0 and 1

# ╔═╡ 67b71def-156a-4ed3-85d3-e7bfd47ffc39
begin
	# Split data, train on 9000 validate on 1000
	tsplit = 9000
	input_data_norm_train = input_data_norm[:, 1:tsplit]
	input_data_norm_valid = input_data_norm[:, tsplit+1:end]
	output_data_norm_train = output_data_norm[:, 1:tsplit]
	output_data_norm_valid = output_data_norm[:, tsplit+1:end]
	nothing
end


# ╔═╡ d562d616-2b98-43fe-b3a3-845177ac3298
begin
	#Loading data into dataloader so we can train on batches at a time
	bs = 512*4
	dataloader = Flux.DataLoader(
	    (input_data = input_data_norm_train,
	     output_data = output_data_norm_train);
	    batchsize = bs
	)
end


# ╔═╡ 38b07258-616b-4ba8-8c2f-dcc2f100a9c5
denormalise(x,min,max) = @. x * (max-min) + min

# ╔═╡ db449505-6d40-4f3b-a9ce-f76ea2a49c06
begin
#After we optimise the parameters, pick the greatest predicted output
idx = argmax(pretrained_model(X_candidate))
best_param_sample = denormalise(X_candidate, inputs_lb, inputs_ub)[:,idx[2]]
inputs_lb .< best_param_sample .< inputs_ub
end

# ╔═╡ f0c3290f-4533-4070-beee-eabfc8abfff8
best_param_sample

# ╔═╡ 9fb68279-d110-4cf6-ac9f-d68972371b4b
begin
	# If you want to see the predicted output for these solutions:
	pred_vals = denormalise(pretrained_model(X_candidate), outputs_lb, outputs_ub)[idx]
end

# ╔═╡ 1cb8c920-4271-4625-abce-76f9baf8be0d
function NNLayer(in_size, out_size, act = identity; bias = false)
    d = Uniform(-1.0 / sqrt(in_size), 1.0 / sqrt(in_size))
    Dense(in_size, out_size, init = (x,y)->rand(d,x,y), act, bias = bias)
end

# ╔═╡ 4605e496-70f9-4e71-8eba-3db00e434145
begin
	act = gelu
	HSIZE = 64
	INSIZE = length(inputs_lb)
	OUTSIZE = length(outputs_lb)
	
	surrogate = Chain(
	    NNLayer(INSIZE, HSIZE, act),
	    SkipConnection(Chain(NNLayer(HSIZE, HSIZE, act),
	                         NNLayer(HSIZE, HSIZE, act)), +),
	    SkipConnection(Chain(NNLayer(HSIZE, HSIZE, act),
	                         NNLayer(HSIZE, HSIZE, act)), +),
	    NNLayer(HSIZE, OUTSIZE)
	)
end


# ╔═╡ 17524f16-3397-4bad-9293-be5a935e73f9
begin
	opt = OptimiserChain(Adam())
	lrs = [1e-3, 5e-4, 3e-4, 1e-4, 5e-5]
	epochs = [100 for i in 1:length(lrs)]
	st = Optimisers.setup(opt, surrogate)
	ps = Flux.params(surrogate) # to use for L2 reg
	lambda_reg = 1e-4
end

# ╔═╡ fc8172cc-9e66-430a-98d0-a4baa37a4c61
begin
loss_t = []
loss_v = []
# Training Loop
for (e, lr) in zip(epochs, lrs)
    for epoch in 1:e
        Optimisers.adjust!(st, lr)
        for batch in dataloader
            x, y = batch.input_data, batch.output_data
            gs = gradient(surrogate) do model
                Flux.mae(model(x), y) + lambda_reg * sum(x_ -> sum(abs2, x_), ps) #l2 reg
            end
            st, surrogate = Optimisers.update(st, surrogate, gs...)
        end
        if epoch % 100 == 0
            surrogate_cpu = surrogate
            l_t = Flux.mae(surrogate_cpu(input_data_norm_train), output_data_norm_train)
            surr_v_pred = denormalise(surrogate_cpu(input_data_norm_valid), outputs_lb, outputs_ub)
            gt_v = denormalise(output_data_norm_valid, outputs_lb, outputs_ub)
            l_v = Flux.mae(surr_v_pred, gt_v)
            push!(loss_t,l_t)
            push!(loss_v,l_v)
            @info "Epoch $epoch lr:$lr Training Loss: $l_t Validation Loss:$l_v"
        end
    end
end
end

# ╔═╡ a2232f58-e65c-4953-a94d-c37d4a27680a
p1 = plot(
	scatter(     
		x = xvals,
		y = loss_t,
		name = "Training Loss"
	),
	Layout(       
		title = "Normalized MAE Loss"
	)
)

# ╔═╡ aaf90e64-3748-4302-a75d-00f8c0a8f067
p2 = plot(
	scatter(
		x = xvals,
		y = loss_v,
		name = "Validation Loss"
	),
	Layout(
		title = "Denormalized MAE Loss",
		xaxis = attr(title = "Epochs")  
	)
)

# ╔═╡ 6fdcaeeb-dc4f-4127-86e1-e54ad2141615

"""
    parallel_coords_plot(params, preds, actual; short_param_labels, colorscale)

Construct a parallel coordinates plot (using PlotlyJS) for `params` (p×N),
`preds` (1×N), and `actual` (1×N). Colors each line by its MAE value.

Parameters
----------
- `params::AbstractMatrix`: A p×N matrix whose columns represent examples.
- `preds::AbstractMatrix`: A 1×N matrix of model predictions (or shape-compatible vector).
- `actual::AbstractMatrix`: A 1×N matrix of ground‐truth values (or shape-compatible vector).
- `short_param_labels::Vector{String}` (optional): Labels for the `p` parameters
  (default: ["param1", "param2", ..., "paramp"]).
- `colorscale`: Plotly colorscale for the lines (default: `[(0,"blue"), (1,"red")]`).

Returns
-------
A PlotlyJS.Plot object, which you can display or save via `PlotlyJS.savefig`.
"""
function parallel_coords_plot(params::AbstractMatrix,
                              preds::AbstractMatrix,
                              actual::AbstractMatrix;
                              short_param_labels::Vector{String} = String[],
                              colorscale = [(0.0, "blue"), (1.0, "red")])

    # 1) Basic checks
    p, N = size(params)
    @assert size(preds, 2) == N "preds must have the same number of columns as params"
    @assert size(actual, 2) == N "actual must have the same number of columns as params"
    @assert size(preds, 1) == 1 "preds should be (1×N) or reshape your data"
    @assert size(actual, 1) == 1 "actual should be (1×N) or reshape your data"

    # If user didn't supply custom labels, create default: "param1", "param2", ...
    if isempty(short_param_labels)
        short_param_labels = ["param$(i)" for i in 1:p]
    else
        @assert length(short_param_labels) == p "short_param_labels must match p (# of rows in params)"
    end

    # 2) Compute MAE for each column
    #    mae_data is a length‐N vector of per‐example MAEs
    mae_data = Float64[]
    for j in 1:N
        push!(mae_data, Flux.mae(preds[:, j], actual[:, j]))
    end

    # 3) Build a DataFrame of all parameters + mae
    df = DataFrame()
    for i in 1:p
        df[!, Symbol("param$(i)")] = params[i, :]  # columns "param1".."paramN"
    end
    df[!, :MAE] = mae_data

    # 4) Construct the "dimensions" array for parcoords
    dimensions_list = Any[]

    #    a) Each param dimension
    for i in 1:p
        param_range = (minimum(df[!, i]), maximum(df[!, i]))
        push!(dimensions_list,
            attr(
                range = param_range,
                label = short_param_labels[i],
                values = df[!, i]
            )
        )
    end

    #    b) The final dimension: MAE
    mae_range = (minimum(mae_data), maximum(mae_data))
    push!(dimensions_list,
        attr(
            range  = mae_range,
            label  = "MAE",
            values = df[!, :MAE]
        )
    )

    # 5) Create the parallel coordinates trace
    mytrace = parcoords(
        line = attr(
            color      = df[!, :MAE],
            colorscale = colorscale
        ),
        dimensions = dimensions_list
    )

    myplot = plot(mytrace)

    return myplot
end


# ╔═╡ 29142a15-5acd-41d0-a273-2437ff805f59
begin
	params_train    = denormalise(input_data_norm_train, inputs_lb, inputs_ub)   # shape (p, N)
	preds_train     = denormalise(surrogate(input_data_norm_train), outputs_lb, outputs_ub)  # shape (1, N)
	actual_train    = output_data[:, 1:9000]  # shape (1, N) grab first 9000
	labels    = ["OroScale", "MtHeight", "MtSize", "MtLon", "MtLat",
	             "TempEqu",  "TempPol",  "TempAtl","TempAzo","ZWind"]
	parallel_coords_plot(params_train, preds_train, actual_train; short_param_labels=labels)
end

# ╔═╡ 0f0a8b66-e312-4f7a-92fd-a1edcf3f343c
begin
	params_valid    = denormalise(input_data_norm_valid, inputs_lb, inputs_ub)   # shape (p, N)
	preds_valid     = denormalise(surrogate(input_data_norm_valid), outputs_lb, outputs_ub)  # shape (1, N)
	actual_valid    = output_data[:, 9000+1:end]  # shape (1, N) grab first 9000
	parallel_coords_plot(params_valid, preds_valid, actual_valid; short_param_labels=labels)
end

# ╔═╡ edc5ea7a-207b-471d-a835-dc6f6db8a63f
begin
	pretrained_preds_train = denormalise(pretrained_model(input_data_norm_train), outputs_lb, outputs_ub) 
	parallel_coords_plot(params_train, pretrained_preds_train, actual_train; short_param_labels=labels)
end

# ╔═╡ 39f1fe53-43a1-4e18-b081-654c1c3a3508
begin
	pretrained_preds_valid = denormalise(pretrained_model(input_data_norm_valid), outputs_lb, outputs_ub) 
	parallel_coords_plot(params_valid, pretrained_preds_valid, actual_valid; short_param_labels=labels)
end

# ╔═╡ e07a3f26-7da0-4dfa-bac3-85a156398e51
function max_precipitation(parameters::AbstractVector)
    parameter_tuple = NamedTuple{PARAMETER_KEYS}(parameters)
    return max_precipitation(parameter_tuple)
end

# ╔═╡ b91ea1a7-b1fa-4988-b8d2-97e9b4b978f9
function max_precipitation(parameters::NamedTuple)

    # define resolution. Use trunc=42, 63, 85, 127, ... for higher resolution, cubically slower
    spectral_grid = SpectralGrid(trunc=31, nlayers=8)

    # Define AquaPlanet ocean, for idealised sea surface temperatures
    # but don't change land-sea mask = retain real ocean basins
    ocean = AquaPlanet(spectral_grid,
                temp_equator=parameters.temperature_equator,
                temp_poles=parameters.temperature_pole)

    initial_conditions = InitialConditions(
        vordiv = ZonalWind(u₀=parameters.zonal_wind),
        temp = JablonowskiTemperature(u₀=parameters.zonal_wind),
        pres = PressureOnOrography(),
        humid = ConstantRelativeHumidity())

    orography = EarthOrography(spectral_grid, scale=parameters.orography_scale)

    # construct model
    model = PrimitiveWetModel(spectral_grid; ocean, initial_conditions, orography)

    # Add rain gauge, locate on Terceira Island
    rain_gauge = RainGauge(spectral_grid, lond=-27.25, latd=38.7)
    add!(model, rain_gauge)

    # Initialize
    simulation = initialize!(model, time=DateTime(2025, 1, 10))

    # Add additional  mountain
    H = parameters.mountain_height
    λ₀, φ₀, σ = parameters.mountain_lon, parameters.mountain_lat, parameters.mountain_size  
    set!(model, orography=(λ,φ) -> H*exp(-spherical_distance((λ,φ), (λ₀,φ₀), radius=360/2π)^2/2σ^2), add=true)

    # set sea surface temperature anomalies
    # 1. Atlantic
    set!(simulation, sea_surface_temperature=
        (λ, φ) -> (30 < φ < 60) && (270 < λ < 360) ? parameters.temperature_atlantic : 0, add=true)

    # 2. Azores
    A = parameters.temperature_azores
    λ_az, φ_az, σ_az = -27.25, 38.7, 4    # location [˚], size [˚] of Azores
    set!(simulation, sea_surface_temperature=
        (λ, φ) -> A*exp(-spherical_distance((λ,φ), (λ_az,φ_az), radius=360/2π)^2/2σ_az^2), add=true)

    # Run simulation for 20 days, maybe longer for more stable statistics? Could be increased to 30, 40, ... days ?
    run!(simulation, period=Day(20))

    # skip first 5 days, as is done in the RainMaker challenge
    RainMaker.skip!(rain_gauge, Day(5))

    # evaluate rain gauge
    lsc = rain_gauge.accumulated_rain_large_scale
    conv = rain_gauge.accumulated_rain_convection
    total_precip = maximum(lsc) + maximum(conv)
    return total_precip
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Flux = "587475ba-b771-5e3f-ad9e-33799f191a9c"
JLSO = "9da8a3cd-07a3-59c0-a743-3fdc52c30d11"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"

[compat]
BenchmarkTools = "~1.5.0"
DataFrames = "~1.7.0"
Distributions = "~0.25.115"
Flux = "~0.16.0"
JLSO = "~2.7.0"
PlutoPlotly = "~0.6.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "7543fd07315f412154f4d5096a016a186819192c"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "e93c42e833e6d4bd28be7b3b56d8deb99fd51f25"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.40"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    DatesExt = "Dates"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "50c3c56a52972d78e8be9fd135bfb91c9574c140"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.1.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgCheck]]
git-tree-sha1 = "680b3b8759bd4c54052ada14e52355ab69e07876"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Atomix]]
deps = ["UnsafeAtomics"]
git-tree-sha1 = "c3b238aa28c1bebd4b5ea4988bebf27e9a01b72b"
uuid = "a9b6321e-bd34-4604-b9c9-b65b8de01458"
version = "1.0.1"

    [deps.Atomix.extensions]
    AtomixCUDAExt = "CUDA"
    AtomixMetalExt = "Metal"
    AtomixoneAPIExt = "oneAPI"

    [deps.Atomix.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    oneAPI = "8f75cd03-7ff8-4ecb-9b8f-daf728133b1b"

[[deps.BSON]]
git-tree-sha1 = "4c3e506685c527ac6a54ccc0c8c76fd6f91b42fb"
uuid = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
version = "0.3.9"

[[deps.BangBang]]
deps = ["Accessors", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires"]
git-tree-sha1 = "e2144b631226d9eeab2d746ca8880b7ccff504ae"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.4.3"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTablesExt = "Tables"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1dff6729bc61f4d49e140da1af55dcd1ac97b2f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.5.0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "SparseInverseSubset", "Statistics", "StructArrays", "SuiteSparse"]
git-tree-sha1 = "4312d7869590fab4a4f789e97bd82f0a04eaaa05"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.72.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "c785dfb1b3bfddd1da557e861b919819b82bbe5b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.27.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7901a6117656e29fa2c74a58adb682f380922c47"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.116"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EnzymeCore]]
git-tree-sha1 = "0cdb7af5c39e92d78a0ee8d0a447d32f7593137e"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.8"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "0a2e5873e9a5f54abb06418d57a8df689336a660"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.2"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "7878ff7172a8e6beedd1dea14bd27c3c6340d361"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.22"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Flux]]
deps = ["Adapt", "ChainRulesCore", "Compat", "EnzymeCore", "Functors", "LinearAlgebra", "MLDataDevices", "MLUtils", "MacroTools", "NNlib", "OneHotArrays", "Optimisers", "Preferences", "ProgressLogging", "Random", "Reexport", "Setfield", "SparseArrays", "SpecialFunctions", "Statistics", "Zygote"]
git-tree-sha1 = "86729467baa309581eb0e648b9ede0aeb40016be"
uuid = "587475ba-b771-5e3f-ad9e-33799f191a9c"
version = "0.16.0"

    [deps.Flux.extensions]
    FluxAMDGPUExt = "AMDGPU"
    FluxCUDAExt = "CUDA"
    FluxCUDAcuDNNExt = ["CUDA", "cuDNN"]
    FluxEnzymeExt = "Enzyme"
    FluxMPIExt = "MPI"
    FluxMPINCCLExt = ["CUDA", "MPI", "NCCL"]

    [deps.Flux.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    NCCL = "3fe64909-d7a1-4096-9b7d-7a0f12cf0f6b"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GPUArrays]]
deps = ["Adapt", "GPUArraysCore", "KernelAbstractions", "LLVM", "LinearAlgebra", "Printf", "Random", "Reexport", "ScopedValues", "Serialization", "Statistics"]
git-tree-sha1 = "86a185085bffbafb11c046f676ee876597b76b66"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "11.2.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.HashArrayMappedTries]]
git-tree-sha1 = "2eaa69a7cab70a52b9687c8bf950a5a93ec895ae"
uuid = "076d061b-32b6-4027-95e0-9a2c6f6d7e74"
version = "0.2.0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "b1c2585431c382e3fe5805874bda6aea90a95de9"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.25"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IRTools]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "950c3717af761bc3ff906c2e8e52bd83390b6ec2"
uuid = "7869d1d1-7146-5819-86e3-90919afe41df"
version = "0.4.14"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JLSO]]
deps = ["BSON", "CodecZlib", "FilePathsBase", "Memento", "Pkg", "Serialization"]
git-tree-sha1 = "7e3821e362ede76f83a39635d177c63595296776"
uuid = "9da8a3cd-07a3-59c0-a743-3fdc52c30d11"
version = "2.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.KernelAbstractions]]
deps = ["Adapt", "Atomix", "InteractiveUtils", "MacroTools", "PrecompileTools", "Requires", "StaticArrays", "UUIDs"]
git-tree-sha1 = "b9a838cd3028785ac23822cded5126b3da394d1a"
uuid = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
version = "0.9.31"
weakdeps = ["EnzymeCore", "LinearAlgebra", "SparseArrays"]

    [deps.KernelAbstractions.extensions]
    EnzymeExt = "EnzymeCore"
    LinearAlgebraExt = "LinearAlgebra"
    SparseArraysExt = "SparseArrays"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Preferences", "Printf", "Unicode"]
git-tree-sha1 = "d422dfd9707bec6617335dc2ea3c5172a87d5908"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "9.1.3"

    [deps.LLVM.extensions]
    BFloat16sExt = "BFloat16s"

    [deps.LLVM.weakdeps]
    BFloat16s = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "TOML"]
git-tree-sha1 = "05a8bd5a42309a9ec82f700876903abce1017dd3"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.34+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MLDataDevices]]
deps = ["Adapt", "Compat", "Functors", "Preferences", "Random"]
git-tree-sha1 = "731113c87c225f372a733522a6d6952ae95ee686"
uuid = "7e8f7934-dd98-4c1a-8fe8-92b47a384d40"
version = "1.6.7"

    [deps.MLDataDevices.extensions]
    MLDataDevicesAMDGPUExt = "AMDGPU"
    MLDataDevicesCUDAExt = "CUDA"
    MLDataDevicesChainRulesCoreExt = "ChainRulesCore"
    MLDataDevicesChainRulesExt = "ChainRules"
    MLDataDevicesComponentArraysExt = "ComponentArrays"
    MLDataDevicesFillArraysExt = "FillArrays"
    MLDataDevicesGPUArraysExt = "GPUArrays"
    MLDataDevicesMLUtilsExt = "MLUtils"
    MLDataDevicesMetalExt = ["GPUArrays", "Metal"]
    MLDataDevicesOneHotArraysExt = "OneHotArrays"
    MLDataDevicesReactantExt = "Reactant"
    MLDataDevicesRecursiveArrayToolsExt = "RecursiveArrayTools"
    MLDataDevicesReverseDiffExt = "ReverseDiff"
    MLDataDevicesSparseArraysExt = "SparseArrays"
    MLDataDevicesTrackerExt = "Tracker"
    MLDataDevicesZygoteExt = "Zygote"
    MLDataDevicescuDNNExt = ["CUDA", "cuDNN"]
    MLDataDevicesoneAPIExt = ["GPUArrays", "oneAPI"]

    [deps.MLDataDevices.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ComponentArrays = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
    FillArrays = "1a297f60-69ca-5386-bcde-b61e274b549b"
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
    MLUtils = "f1d291b0-491e-4a28-83b9-f70985020b54"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    OneHotArrays = "0b1bfda6-eb8a-41d2-88d8-f5af5cad476f"
    Reactant = "3c362404-f566-11ee-1572-e11a4b42c853"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"
    oneAPI = "8f75cd03-7ff8-4ecb-9b8f-daf728133b1b"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MLUtils]]
deps = ["ChainRulesCore", "Compat", "DataAPI", "DelimitedFiles", "FLoops", "NNlib", "Random", "ShowCases", "SimpleTraits", "Statistics", "StatsBase", "Tables", "Transducers"]
git-tree-sha1 = "7940c0af802586b97009f254aa6065000a16fa1d"
uuid = "f1d291b0-491e-4a28-83b9-f70985020b54"
version = "0.4.5"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Memento]]
deps = ["Dates", "Distributed", "Requires", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "bb2e8f4d9f400f6e90d57b34860f6abdc51398e5"
uuid = "f28f55f0-a522-5efc-85c2-fe41dfb9b2d9"
version = "1.4.1"

[[deps.MicroCollections]]
deps = ["Accessors", "BangBang", "InitialValues"]
git-tree-sha1 = "44d32db644e84c75dab479f1bc15ee76a1a3618f"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.2.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NNlib]]
deps = ["Adapt", "Atomix", "ChainRulesCore", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "bdc9d30f151590aca0af22690f5ab7dc18a551cb"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.9.27"

    [deps.NNlib.extensions]
    NNlibAMDGPUExt = "AMDGPU"
    NNlibCUDACUDNNExt = ["CUDA", "cuDNN"]
    NNlibCUDAExt = "CUDA"
    NNlibEnzymeCoreExt = "EnzymeCore"
    NNlibFFTWExt = "FFTW"
    NNlibForwardDiffExt = "ForwardDiff"

    [deps.NNlib.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    cuDNN = "02a925ec-e4fe-4b08-9a7e-0d78e3d38ccd"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "030ea22804ef91648f29b7ad3fc15fa49d0e6e71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.3"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OneHotArrays]]
deps = ["Adapt", "ChainRulesCore", "Compat", "GPUArraysCore", "LinearAlgebra", "NNlib"]
git-tree-sha1 = "c8c7f6bfabe581dc40b580313a75f1ecce087e27"
uuid = "0b1bfda6-eb8a-41d2-88d8-f5af5cad476f"
version = "0.2.6"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optimisers]]
deps = ["ChainRulesCore", "Functors", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "53ff746a3a2b232a37dbcd262ac8bbb2b18202b8"
uuid = "3bd65402-5787-11e9-1adc-39752487f4e2"
version = "0.4.4"
weakdeps = ["Adapt", "EnzymeCore"]

    [deps.Optimisers.extensions]
    OptimisersAdaptExt = ["Adapt"]
    OptimisersEnzymeCoreExt = "EnzymeCore"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Artifacts", "ColorSchemes", "Colors", "Dates", "Downloads", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "Pkg", "PlotlyBase", "PrecompileTools", "Reexport", "ScopedValues", "Scratch", "TOML"]
git-tree-sha1 = "9ebe25fc4703d4112cc418834d5e4c9a4b29087d"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.6.2"

    [deps.PlutoPlotly.extensions]
    PlotlyKaleidoExt = "PlotlyKaleido"
    UnitfulExt = "Unitful"

    [deps.PlutoPlotly.weakdeps]
    PlotlyKaleido = "f2990250-8cf9-495f-b13a-cce12b45703c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.PtrArrays]]
git-tree-sha1 = "77a42d78b6a92df47ab37e177b2deac405e1c88f"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "cda3b045cf9ef07a08ad46731f5a3165e56cf3da"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.1"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.ScopedValues]]
deps = ["HashArrayMappedTries", "Logging"]
git-tree-sha1 = "1147f140b4c8ddab224c94efa9569fc23d63ab44"
uuid = "7e506255-f358-4e82-b7e4-beb19740aa63"
version = "1.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShowCases]]
git-tree-sha1 = "7f534ad62ab2bd48591bdeac81994ea8c445e4a5"
uuid = "605ecd9f-84a6-4c9e-81e2-4798472b76a3"
version = "0.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SparseInverseSubset]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "52962839426b75b3021296f7df242e40ecfc0852"
uuid = "dc90abb0-5640-4711-901d-7e5b23a2fada"
version = "0.1.2"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "47091a0340a675c738b1304b58161f3b0839d454"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.10"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "9537ef82c42cdd8c5d443cbc359110cbb36bae10"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.21"
weakdeps = ["Adapt", "GPUArraysCore", "KernelAbstractions", "LinearAlgebra", "SparseArrays", "StaticArrays"]

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Transducers]]
deps = ["Accessors", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "SplittablesBase", "Tables"]
git-tree-sha1 = "7deeab4ff96b85c5f72c824cae53a1398da3d1cb"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.84"

    [deps.Transducers.extensions]
    TransducersAdaptExt = "Adapt"
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnsafeAtomics]]
git-tree-sha1 = "b13c4edda90890e5b04ba24e20a310fbe6f249ff"
uuid = "013be700-e6cd-48c3-b4a1-df204f14c38f"
version = "0.3.0"
weakdeps = ["LLVM"]

    [deps.UnsafeAtomics.extensions]
    UnsafeAtomicsLLVM = ["LLVM"]

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zygote]]
deps = ["AbstractFFTs", "ChainRules", "ChainRulesCore", "DiffRules", "Distributed", "FillArrays", "ForwardDiff", "GPUArrays", "GPUArraysCore", "IRTools", "InteractiveUtils", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NaNMath", "PrecompileTools", "Random", "Requires", "SparseArrays", "SpecialFunctions", "Statistics", "ZygoteRules"]
git-tree-sha1 = "0b3c944f5d2d8b466c5d20a84c229c17c528f49e"
uuid = "e88e6eb3-aa80-5325-afca-941959d7151f"
version = "0.6.75"

    [deps.Zygote.extensions]
    ZygoteColorsExt = "Colors"
    ZygoteDistancesExt = "Distances"
    ZygoteTrackerExt = "Tracker"

    [deps.Zygote.weakdeps]
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    Distances = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "27798139afc0a2afa7b1824c206d5e87ea587a00"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─a9c6cecd-0013-42b5-b7cc-0362bdb4a663
# ╟─7e34c342-91dd-4c15-8a4c-0d6a61093371
# ╠═97a4183e-f670-47b6-a318-4254f34cac32
# ╟─62dba600-c90b-11ef-009f-07f05e9e18cc
# ╠═a26dfb20-2644-4657-b5b5-864a32ac3053
# ╟─af87db2b-31af-4aa9-b1a9-9d335bcde0e8
# ╟─0c685bdf-edaf-4116-8f05-7e11383e3200
# ╠═08523ab9-5c66-4ceb-9c12-96ad1e5e9c01
# ╠═c13793fa-30ea-43e7-a416-7651877c140b
# ╠═f1ee0a37-cae6-4223-a5b0-aea809fca768
# ╟─e57ee8db-fa0f-422f-8c0a-7ee612cc732a
# ╟─106e4cd3-e376-4c5b-b0eb-90952d79c760
# ╠═77e3c694-1731-4e14-98e8-b09cfbaa9f94
# ╠═2c0b78e5-817f-4eea-b0e1-41d5cfcc8cc3
# ╠═f52bba40-fd7d-4646-82d3-cc6420219916
# ╠═73ae3f0f-7bbf-44d5-8522-acce0c9afc3a
# ╠═67b71def-156a-4ed3-85d3-e7bfd47ffc39
# ╠═d562d616-2b98-43fe-b3a3-845177ac3298
# ╟─b5b681a6-f020-41c3-a690-51a9f730515b
# ╠═4605e496-70f9-4e71-8eba-3db00e434145
# ╟─31e025b7-8b08-4adb-bc55-c4ebd72a931b
# ╠═17524f16-3397-4bad-9293-be5a935e73f9
# ╠═22dad96e-c4c8-4ff8-a4d7-f83accb02b72
# ╠═fc8172cc-9e66-430a-98d0-a4baa37a4c61
# ╠═70bfa951-27b3-4688-83c2-dd5be03ff244
# ╠═a2232f58-e65c-4953-a94d-c37d4a27680a
# ╠═aaf90e64-3748-4302-a75d-00f8c0a8f067
# ╠═607321b9-fd3f-41bf-91d0-2e7d099afc6c
# ╠═29142a15-5acd-41d0-a273-2437ff805f59
# ╟─df9b4742-f68e-4ff3-b501-ecbe3199d300
# ╠═0f0a8b66-e312-4f7a-92fd-a1edcf3f343c
# ╟─ca4456ab-dc1c-402c-98e1-07b7665bf86d
# ╠═10e14cb4-ab21-4eda-b106-bc785c4ea859
# ╟─8b7d5e71-d10b-493e-addf-6df24971e385
# ╠═edc5ea7a-207b-471d-a835-dc6f6db8a63f
# ╟─60770e3b-83a4-4a8a-a631-b92b0e0242d7
# ╠═39f1fe53-43a1-4e18-b081-654c1c3a3508
# ╟─7c27f6c2-4869-4cee-9f49-b8530dde3c6f
# ╠═f0441313-65d5-4a75-8499-2f17db7668a7
# ╠═796e1216-fc68-41f0-b47b-1d638dff6bf4
# ╠═db449505-6d40-4f3b-a9ce-f76ea2a49c06
# ╠═f0c3290f-4533-4070-beee-eabfc8abfff8
# ╠═9fb68279-d110-4cf6-ac9f-d68972371b4b
# ╟─199c8bf6-4ef7-43ee-8a44-d168dbf94304
# ╟─1f154133-62d4-47f6-8232-b5a35f927378
# ╟─cae4b323-475f-44f4-8586-d3b41e1ae8b1
# ╟─29dd72ce-1145-436d-892c-c3b5d8693d69
# ╟─5b4c9729-aafe-4231-9429-a11ada56baf1
# ╟─18aeca69-12d4-4414-932e-f87ba42f79a6
# ╟─38b07258-616b-4ba8-8c2f-dcc2f100a9c5
# ╟─1cb8c920-4271-4625-abce-76f9baf8be0d
# ╟─6fdcaeeb-dc4f-4127-86e1-e54ad2141615
# ╟─e07a3f26-7da0-4dfa-bac3-85a156398e51
# ╟─b91ea1a7-b1fa-4988-b8d2-97e9b4b978f9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
