### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 71fcff82-cb70-11ef-18ff-3d0440ed8e1c
begin
	using Pkg
	Pkg.add("BinStatistics")
	Pkg.add("DataFrames")
	Pkg.add("Statistics")
	Pkg.add("CairoMakie")
	Pkg.add("Pluto")
	using BinStatistics
	using DataFrames
	using Statistics
	using CairoMakie
end

# ╔═╡ 1167ad2f-ee64-465b-b8bc-53dea3ca5641
md"
# BinStatistics.jl Demo
"

# ╔═╡ ecdde429-e34c-4f06-8417-a1ea4cf4c924
md"""

##### BinStatistics.jl exports a single *binstats* function:

    binstats(
		df, 
		axis_col, 
		axis_edges, 
		bin_col; 
        grp_function = nrow, 
		col_function = mean, 
		missing_bin = false
	)
    
Returns a DataFrame containing function values for binned variables of *df*.

#### Arguments
- *axis_col*: binning axes column(s)
- *axis_edges*: bin edges for *axis_col*
- *bin_col*: column variable(s) to be binned
- *grp_function* = nrow: column independent funciton(s) to be applied at group level
- *var_function* = mean: column dependent funciton(s) to be applied to *bin_col* at group level
- *missing_bins = false*: include missing bins
"""

# ╔═╡ 565fc12c-7b6b-4b29-bd89-0fd69a586994
md"""
## Load packages
"""

# ╔═╡ a49c141a-418d-4baf-80c7-be19f0e2420b
md"""
## Create synthetic data
"""

# ╔═╡ 3127244c-14c3-4864-a7f0-fc24a9b5a5c0
begin
    n = 100000
    df = DataFrame()
    df[!, :x] = rand(n) .* 20
    df[!, :y] = rand(n) .* 20
    df[!, :v1] = cos.(df.x) .+ randn(n) * 3
    df[!, :v2] = cos.(df.x .- df.y) .+ sin.(df.x .+ df.y) .+ randn(n) * 3
    df[!, :v3] = df.v1 .+ df.v2
end;

# ╔═╡ a60499ce-3215-453e-a323-b14c670a8831
md"## Example 1: calculate count/nrow and mean of *v1* binned according to *x*"

# ╔═╡ 350bc0a9-2c1a-4b20-820d-175ecd7022c8
begin
	axis_var1 = :x
	axis_edges1 = 0:0.1:20
	bin_var1 = :v1
	df1 = binstats(df, axis_var1, axis_edges1, bin_var1)
end

# ╔═╡ 3675dad9-975c-48f5-9ea3-ca1213464fc6
begin
    fig1 = Figure()
    Axis(fig1[1, 1], title="raw data")
    scatter!(fig1[1, 1], df.x, df.v1)
    Axis(fig1[1, 2], title="binned data")
    scatter!(fig1[1, 2], df1[:, 1], df1.v1_mean)
    fig1
end

# ╔═╡ a7dcc2d3-bf35-4658-afe8-3b3de0e88d3e
md"## Example 2: calculate count/nrow and mean of *v1* and *v2* binned according to *x*"

# ╔═╡ a32b169b-3b2a-4a3c-b337-41f56f37faac
begin
	axis_var2 = :x
	axis_edges2 = 0:0.1:20
	bin_var2 = [:v1, :v2]
	
	df2 = binstats(df, axis_var2 , axis_edges2, bin_var2)
end


# ╔═╡ ec6ce34d-7dda-4f94-904d-c0f837a32740
begin fig2 = Figure()
    Axis(fig2[1, 1], title="raw data")
    scatter!(fig2[1, 1], df.x, df.v1)
    scatter!(fig2[1, 1], df.x, df.v2)
    Axis(fig2[1, 2], title="binned data")
    scatter!(fig2[1, 2], df2[:, 1], df2.v1_mean, label="v1")
    scatter!(fig2[1, 2], df2[:, 1], df2.v2_mean, label="v2")
    axislegend()
    fig2
end

# ╔═╡ 286e3ae9-6315-4e63-8416-41f2e6098c09
md"## Example 3: calculate count/nrow, mean, medain and std of *v1* binned according to *x"

# ╔═╡ fb078840-36fc-4a32-af9d-36b42ce5d6d1
begin
	axis_var3 = :x
	axis_edges3 = 0:0.1:20
	bin_var3 = :v1
	col_function3 = [mean, median, std]
	
	df3 = binstats(df, axis_var3, axis_edges3, bin_var3; col_function = col_function3)
end

# ╔═╡ 320cbcec-af9a-4cc4-909b-13f5d620410c
begin
    fig3 = Figure()
    Axis(fig3[1, 1], title="raw data")
    scatter!(fig3[1, 1], df.x, df.v1)
    Axis(fig3[1, 2], title="binned data")
    scatter!(fig3[1, 2], df3[:, 1], df3.v1_mean, label="mean")
    scatter!(fig3[1, 2], df3[:, 1], df3.v1_median, label="median")
    scatter!(fig3[1, 2], df3[:, 1], df3.v1_std, label="std")
    axislegend()
    fig3
end

# ╔═╡ 8ec188ce-5afe-4a26-b443-865d246d734e
md"## Example 4: calculate count/nrow and mean of *v2* binned according to *y* and *x*"

# ╔═╡ dc9ebb91-9c98-426e-8d3a-89d5c3f6035c
begin
	axis_var4 = [:y, :x]
	axis_edges4 = [0:0.2:20, 0:0.2:20]
	bin_var4 = :v2
	
	df4 = binstats(df, axis_var4, axis_edges4, bin_var4; missing_bins = true)
end

# ╔═╡ 661923eb-8483-4113-ae73-a0e5ea884648
begin
    fig4 = Figure()
    Axis(fig4[1, 1], title="raw data")
    scatter!(fig4[1, 1], df.y, df.x, color=df.v2, colormap=:thermal, markersize=1)
    xlims!(0, 20)
    ylims!(0, 20)
    Axis(fig4[1, 2], title="binned data")
    heatmap!(fig4[1, 2], unique(df4[:, 1]), unique(df4[:, 2]),
        reshape(df4.v2_mean, length(unique(df4[:, 2])), length(unique(df4[:, 1]))),
        colormap=:thermal)
    fig4
end

# ╔═╡ 7066d731-2990-4020-99dc-960880cdbbd5
md"## Example 5: calculate median of *v2* binned according to *y* and *x* using non-uniform axis_edges"


# ╔═╡ b9857e5d-e345-464a-ad4c-e55621f15b43
begin
	axis_var5 = [:y, :x]
	axis_edges5 = [(0:0.5:4.5) .^ 2, (0:0.5:4.5) .^ 2]
	bin_var5 = :v2
	missing_bins5 = true
	col_function5 = median
	
	df5 = binstats(df, axis_var5, axis_edges5, bin_var5; col_function = col_function5, missing_bins = true)
end


# ╔═╡ 82d63e53-f42f-41a9-a883-de366bbfff2d

begin
    fig5 = Figure()
    Axis(fig5[1, 1], title="raw data")
    scatter!(fig5[1, 1], df.y, df.x, color=df.v2, colormap=:thermal, markersize=1)
    xlims!(0, 20)
    ylims!(0, 20)
    Axis(fig5[1, 2], title="binned data")
    heatmap!(fig5[1, 2], unique(df5[:, 1]), unique(df5[:, 2]),
        reshape(df5.v2_median, length(unique(df5[:, 2])), length(unique(df5[:, 1]))), colormap=:thermal)
    fig5
end

# ╔═╡ 05bb54db-d381-49d6-8372-be9a51ca0950
md"## Example 6: apply custom function to *v2*, binned according to *y* and *x* create a median absolute deviation function"

# ╔═╡ a4b693d0-a4cc-46fa-ac61-2c642cf2cca5
begin
	function mad(x)
	    median(abs.(x .- median(x))) 
	end
	# binstats also accepts anonymous functions but the output will be assinged a generic name
	
	axis_var6 = [:y, :x]
	axis_edges6 = [0:1:20, 0:1:20]
	bin_var6 = :v2
	missing_bins6 = true
	col_function6 = mad
	missing_bins = true

	# apply to grouped data
	df6 = binstats(df, axis_var6, axis_edges6, bin_var6; col_function = col_function6, missing_bins = true)
end



# ╔═╡ 3cac9e3c-da45-49c8-a523-9b8b67e4836b
begin
    fig6 = Figure()
    Axis(fig6[1, 1], title="raw data")
    scatter!(fig6[1, 1], df.y, df.x, color=df.v2, colormap=:thermal, markersize=1)
    xlims!(0, 20)
    ylims!(0, 20)
    Axis(fig6[1, 2], title="binned data")
    heatmap!(fig6[1, 2], unique(df6[:, 1]), unique(df6[:, 2]),
        reshape(df6.v2_mad, length(unique(df6[:, 2])), length(unique(df6[:, 1]))),
        colormap=:thermal)
    fig6
end

# ╔═╡ Cell order:
# ╟─1167ad2f-ee64-465b-b8bc-53dea3ca5641
# ╟─ecdde429-e34c-4f06-8417-a1ea4cf4c924
# ╟─565fc12c-7b6b-4b29-bd89-0fd69a586994
# ╠═71fcff82-cb70-11ef-18ff-3d0440ed8e1c
# ╟─a49c141a-418d-4baf-80c7-be19f0e2420b
# ╠═3127244c-14c3-4864-a7f0-fc24a9b5a5c0
# ╠═a60499ce-3215-453e-a323-b14c670a8831
# ╠═350bc0a9-2c1a-4b20-820d-175ecd7022c8
# ╟─3675dad9-975c-48f5-9ea3-ca1213464fc6
# ╟─a7dcc2d3-bf35-4658-afe8-3b3de0e88d3e
# ╠═a32b169b-3b2a-4a3c-b337-41f56f37faac
# ╟─ec6ce34d-7dda-4f94-904d-c0f837a32740
# ╟─286e3ae9-6315-4e63-8416-41f2e6098c09
# ╠═fb078840-36fc-4a32-af9d-36b42ce5d6d1
# ╟─320cbcec-af9a-4cc4-909b-13f5d620410c
# ╟─8ec188ce-5afe-4a26-b443-865d246d734e
# ╠═dc9ebb91-9c98-426e-8d3a-89d5c3f6035c
# ╟─661923eb-8483-4113-ae73-a0e5ea884648
# ╟─7066d731-2990-4020-99dc-960880cdbbd5
# ╠═b9857e5d-e345-464a-ad4c-e55621f15b43
# ╟─82d63e53-f42f-41a9-a883-de366bbfff2d
# ╟─05bb54db-d381-49d6-8372-be9a51ca0950
# ╠═a4b693d0-a4cc-46fa-ac61-2c642cf2cca5
# ╟─3cac9e3c-da45-49c8-a523-9b8b67e4836b
