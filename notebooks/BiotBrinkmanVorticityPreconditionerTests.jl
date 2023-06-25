### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 38e1939a-3c91-4b7b-88fd-ae7e01eda602
# Run notebook in Pluto.jl "backward compatibility mode" 
# using the so-called "global environment" pattern
# see https://github.com/fonsp/Pluto.jl/wiki/%F0%9F%8E%81-Package-management 
# for more details
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
end

# ╔═╡ 299c6e30-9107-4039-a11c-32ed6d9b1460
using DrWatson

# ╔═╡ 52a4a234-a99c-4e0d-b1e7-f54ddc11ff45
using DataFrames

# ╔═╡ 08a6ab1f-ef85-41a3-bb1a-933d578e4f93
using BSON

# ╔═╡ 696ffce0-e025-409c-8628-48b06fc7fd38
using PlutoUI

# ╔═╡ e82ad202-94eb-4f79-8aa3-7dd40375328b
using Plots

# ╔═╡ 63536890-4206-48eb-897e-e868e371c721
using LaTeXStrings

# ╔═╡ 106ded38-479d-4961-a127-24c2b416aa9d
md"""## BiotBrinkmanVorticityPreconditionerTests
"""

# ╔═╡ 4b25a735-1391-4985-a761-6764b20985c5
gr()

# ╔═╡ 69205551-a348-4d20-8b1e-0860f6c7eafe
experiment_data_directory="BiotBrinkmanVorticityPreconditionerTests";

# ╔═╡ efe130ca-eb74-4cb4-af77-f6be72f048f6
commit_dirs=readdir(datadir(experiment_data_directory)); 

# ╔═╡ 2795d244-a4ce-4e0c-96f4-469a432b96f6
md"""
Select commit ID data directory: $(@bind commitID Select(commit_dirs)) 
"""

# ╔═╡ 7436ace4-c30a-474d-b0cd-195f7fd44d77
df = collect_results(datadir(experiment_data_directory,commitID));

# ╔═╡ 594ac663-8106-4bec-a5c4-429c06795044
rename!(df,:prec_variant => :B);

# ╔═╡ 602c6a8b-5a9d-4aba-9a73-c7af073feed0
cols_to_filter = [:μ,:λ,:ν,:κ,:α,:c_0,:ndofs,:niters,:B];

# ╔═╡ 82f442f6-9b3b-46ce-9c91-bb3e54abd526
df_filtered = df[:,cols_to_filter];

# ╔═╡ e4e03f80-42b4-459a-b75a-1e1364b3b659
df_filtered_cols = Dict(pairs(eachcol(df_filtered)));

# ╔═╡ 414e8521-662f-40f9-9c9e-def65d438fc7
params_possible_values=
	    Dict([k=>unique(df_filtered_cols[k]) 
			    for k in keys(df_filtered_cols)]);

# ╔═╡ 17e099ff-caec-4d52-8de8-c3d9f110cf0c


# ╔═╡ e6574182-c313-4f09-b9a7-ad09751ad5ef
md"""
Select the parameter-value combination that you want to visualize!:

μ:   $(@bind μval Select(params_possible_values[:μ]))
λ:   $(@bind λval Select(params_possible_values[:λ]))
ν:   $(@bind νval Select(params_possible_values[:ν]))
κ:   $(@bind κval Select(params_possible_values[:κ]))
α:   $(@bind αval Select(params_possible_values[:α])) 
c0: $(@bind c0val Select(params_possible_values[:c_0]))

Customize visualization

legend position: $(@bind lposition Select([:right, :left, :top, :bottom, :inside, :best, :legend, :topright, :topleft, :bottomleft, :bottomright, :outertopleft];default=:topleft))

autoxlims: $(@bind autoxlims CheckBox(;default=true)) 
xlimleft: $(@bind xliml TextField((2,1);default="0.0"))
xlimright: $(@bind xlimr TextField((2,1);default="1.0"))

autoylims: $(@bind autoylims CheckBox(;default=false))
ylimbottom: $(@bind ylimb TextField((2,1);default="0.0"))
ylimtop: $(@bind ylimt TextField((2,1);default="500.0"))

logx: $(@bind logxval CheckBox(;default=true)) 
logy: $(@bind logyval CheckBox())

"""

# ╔═╡ 5748a503-4534-462a-b485-386d68fe857d
function generate_labels(params_possible_values,keys_to_filter)
  dl=dict_list(params_possible_values)
  labels=Vector{String}(undef,length(dl))
  for (i,d) in enumerate(dl)
	label=""
    if (haskey(d,:c_0) && haskey(d,:κ))
	  if (d[:c_0]==1.0 && d[:κ]==1.0)
		  label=L"c_0=1 \ \kappa=1"
	  elseif (d[:c_0]==1.0 && d[:κ]==1.0e-08)
		  label=L"c_0=1 \ \kappa=10^{-8}"
	  elseif (d[:c_0]==1.0e-08 && d[:κ]==1.0)
		  label=L"c_0=10^{-8} \ \kappa=1"
	  elseif (d[:c_0]==1.0e-08 && d[:κ]==1.0e-08)
		  label=L"c_0=10^{-8} \ \kappa=10^{-8}"	  
	  end 
	else 	  
		for (key,val) in d
		  if !(key in keys_to_filter)
			label=label * "$(key)=$(val)"	
		  end 	  
		end
	end
	labels[i]=label  
  end
  labels
end 

# ╔═╡ 2c7a1e9b-ebe7-4631-a14d-0c2ea0d7293b
function get_x_y(xparam, yparam, ffilter, df)
  df_filtered = filter(ffilter,df)
  df_filtered_cols = Dict(pairs(eachcol(df_filtered)))
  @assert xparam in keys(df_filtered_cols)
  @assert yparam in keys(df_filtered_cols)
	
  params_possible_values=
	    Dict([k=>unique(df_filtered_cols[k]) 
			    for k in keys(df_filtered_cols) if k != xparam && k != yparam])
	
  dl=dict_list(params_possible_values)

  # The following code is general enough so that for 
  # fixed (xparam, yparam) there might be several 
  # possible combinations for the rest of parameters
  # after applying ffilter. In such a case we generate
  # as many curves as combinations of the rest of 
  # parameter values.
  xall=[]
  yall=[]	
  for d in dl
      function f(a...)
		  equal=all(a .== values(d)) 
	  end 
	  ffilter_current_d = collect(keys(d))=>f
	  df_tmp=filter(ffilter_current_d,df_filtered)
	  sort!(df_tmp,[xparam,])
      x = df_tmp[!,xparam]
      y = df_tmp[!,yparam]
      push!(xall,x)
	  push!(yall,y)
  end
  (xall,yall,params_possible_values)
end

# ╔═╡ 9478d88e-bc39-40b2-be09-3789e4fff51e
function plot_xparam_versus_yparam(xparam,yparam,xaxis,yaxis,ffilter,df;
                                   autoxlims=true,
                                   autoylims=true,
                                   title="",
                                   titlefontsize=14,
                                   legendfontsize=10,
                                   xliml=0.0,
                                   xlimr=1.0,
                                   ylimb=0.0,
                                   ylimt=500.0,
                                   markersize=4,
                                   extra_keys_to_filter=Symbol[],
                                   ylabel="$xparam",
                                   xlabel="$yparam",
                                   xtickfontsize=5,
                                   ytickfontsize=5, 
                                   xlabelfontsize=5,
                                   ylabelfontsize=5)
  f = plot()
  x,y,params_possible_values = get_x_y(xparam,yparam,ffilter,df)
  keys_to_filter=deepcopy(first(ffilter))
  push!(keys_to_filter, extra_keys_to_filter...)	
  labels=generate_labels(params_possible_values,keys_to_filter)
  @assert length(x)==length(y)
  @assert length(labels)==length(y)	
  for (xi,yi,li) in zip(x,y,labels)
    plot!(xi,yi,xaxis=xaxis,yaxis=yaxis,
		  title=title,titlefontsize=titlefontsize,
		  legendfontsize=legendfontsize,
		  label=li,markershape=:auto,markersize=markersize,
	      xlabel=xlabel,ylabel=ylabel,xtickfontsize=xtickfontsize,
	      ytickfontsize=ytickfontsize,xlabelfontsize=xlabelfontsize,
	      ylabelfontsize=ylabelfontsize,thickness_scaling=1,
		  legendmarkersize=8)
  end 
  plot!(xlabel=xlabel,ylabel=ylabel,legend=lposition)
  if (!autoxlims)
    xlims!((xliml,xlimr))
  end	  
  if (!autoylims)
    ylims!((ylimb,ylimt))
  end	
  f
end

# ╔═╡ 9880a404-126f-451c-99a9-4787ac816ec5
function generate_mxn_grid_plot(grid_param,
	                            layout,
	                            xparam,
	                            yparam,
	                            ffilter,
	                            df;
	                            title="",
								titlefontsize=titlefontsize,
							    legendfontsize=legendfontsize,
								xtickfontsize=5,
								ytickfontsize=5,
                                size=size,
                                markersize=markersize,
                                xlabel="",
                                ylabel="",
                                xlabelfontsize=5,
			                    ylabelfontsize=5)
	@assert !(grid_param in first(ffilter))
	_,_,params_possible_values = get_x_y(xparam,yparam,ffilter,df)
	@assert grid_param in keys(params_possible_values)
	grid_param_values=params_possible_values[grid_param]
    plots=Vector{Any}(undef,length(grid_param_values))
	for (i,val) in enumerate(grid_param_values)
		ffilter_grid_param=[grid_param]=>(grid_param)->(grid_param==val);
		df_filtered_grid_param = filter(ffilter_grid_param,df)
		println(title)
		if (grid_param==:B && val==:B1)
		   title_subplot=L"\mathcal{B}=\mathcal{B}1 \ " * title
		elseif (grid_param==:B && val==:B4)
		   title_subplot=L"\mathcal{B}=\mathcal{B}3 \ " * title
		end 
		plots[i]=
	        plot_xparam_versus_yparam(xparam,yparam,
		                  (logxval ? :log10 : :none),
		                  (logyval ? :log10 : :none),
	                      ffilter,
		                  df_filtered_grid_param;
	                      autoxlims=autoxlims,
				          autoylims=autoylims,
						  title=title_subplot,
						  titlefontsize=titlefontsize,
						  legendfontsize=legendfontsize,
						  markersize=markersize,
	                      xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
						  ylimb=parse(Float64, ylimb),ylimt=parse(Float64, ylimt),
			              extra_keys_to_filter=[grid_param],
			              xlabel=xlabel,
			              ylabel=ylabel,
			              xtickfontsize=xtickfontsize,
			              ytickfontsize=ytickfontsize,
				          xlabelfontsize=xlabelfontsize,
			              ylabelfontsize=ylabelfontsize,
			              )
	end 
	plot(plots..., layout = layout, size=size)
end 

# ╔═╡ af9c277b-6ad2-4be0-abb7-5ad505aa7649
begin
	cols_to_filter_plot1 = [:μ,:λ,:ν,:κ,:α,:c_0,:ndofs,:niters,:B];
	if (μval==1.0)
		μvals="1"
	elseif (μval==1.0e-08)
		μvals="10^{-8}"
	end
	if (λval==1.0)
		λvals="1"
	elseif (λval==1.0e08)
		λvals="10^{8}"
	end
	if (αval==1.0)
		αvals="1"
	elseif (αval==1.0e-08)
		αvals="10^{-8}"
	end
	title_plot1=L"\ \ \mu=%$(μvals) \ \ \lambda=%$(λvals) \ \ \alpha=%$(αvals) \ \ \nu=1"
	ffilter_plot1=[:μ,:λ,:α,:ν]=>(μ,λ,α,ν)->(μ==μval && λ==λval && α==αval && ν==1);
	df_filtered_plot1 = df[:,cols_to_filter_plot1];
	layout = @layout [a b];
	plt=generate_mxn_grid_plot(:B, layout,
		                   :ndofs, :niters,
		                   ffilter_plot1, df_filtered_plot1;
	                       title=title_plot1,
	                       titlefontsize=12,
	                       legendfontsize=8,
						   xtickfontsize=12,
						   ytickfontsize=12,
						   xlabelfontsize=16,
			               ylabelfontsize=16,
						   markersize=6,
	                       size=(800,500),
	                       xlabel=L"\mathrm{DoF}",
						   ylabel=L"\mathrm{\#iterations}")
    savefig(plt,"iters_versus_dofs_lambda_1_v1e-08.pdf")
	plt
end

# ╔═╡ Cell order:
# ╟─106ded38-479d-4961-a127-24c2b416aa9d
# ╟─38e1939a-3c91-4b7b-88fd-ae7e01eda602
# ╠═299c6e30-9107-4039-a11c-32ed6d9b1460
# ╠═52a4a234-a99c-4e0d-b1e7-f54ddc11ff45
# ╠═08a6ab1f-ef85-41a3-bb1a-933d578e4f93
# ╠═696ffce0-e025-409c-8628-48b06fc7fd38
# ╠═e82ad202-94eb-4f79-8aa3-7dd40375328b
# ╠═4b25a735-1391-4985-a761-6764b20985c5
# ╠═63536890-4206-48eb-897e-e868e371c721
# ╠═69205551-a348-4d20-8b1e-0860f6c7eafe
# ╠═efe130ca-eb74-4cb4-af77-f6be72f048f6
# ╠═2795d244-a4ce-4e0c-96f4-469a432b96f6
# ╠═7436ace4-c30a-474d-b0cd-195f7fd44d77
# ╠═594ac663-8106-4bec-a5c4-429c06795044
# ╠═602c6a8b-5a9d-4aba-9a73-c7af073feed0
# ╠═82f442f6-9b3b-46ce-9c91-bb3e54abd526
# ╠═e4e03f80-42b4-459a-b75a-1e1364b3b659
# ╠═414e8521-662f-40f9-9c9e-def65d438fc7
# ╟─17e099ff-caec-4d52-8de8-c3d9f110cf0c
# ╟─e6574182-c313-4f09-b9a7-ad09751ad5ef
# ╠═af9c277b-6ad2-4be0-abb7-5ad505aa7649
# ╠═9880a404-126f-451c-99a9-4787ac816ec5
# ╠═9478d88e-bc39-40b2-be09-3789e4fff51e
# ╠═5748a503-4534-462a-b485-386d68fe857d
# ╟─2c7a1e9b-ebe7-4631-a14d-0c2ea0d7293b
