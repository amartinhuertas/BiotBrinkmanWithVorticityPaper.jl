
using DrWatson

# The following macro call let us execute the script with
# the project's environment even if we ran the julia REPL
# without the --project=... flag
@quickactivate "BiotBrinkmanWithVorticityPaper"

function replace_strings_by_symbols(r)
  d = Dict{Symbol,Any}()
  for (k,v) in r
    d[Symbol(k)] = v
  end
  d
end

# Define parameter-value combinations for the experiment.
# Parameter-value combinations s.t. corresponding results
# are already available in the data folder are not re-run.
# You have to eliminate
# them to be re-run

function generate_param_dicts()
   params = Dict{Symbol,Any}(
     :nk           => [1, 2, 3, 4, 5],
     :order        => [0],
     :μ            => [1.,1.e8],
     :λ            => [1.,1.e8],
     :ν            => [1.,1.e-8],
     :κ            => [1.,1.e-8],
     :α            => [1.],
     :c_0          => [1.,1.e-8],
     :prec_variant => [:B1,:B2,:B3],
   )
   dict_list(params)
end

# Defines the Driver module with the driver(...) function inside
# The computational heavy stuff and the actual code of the experiment
# at hand is here.
include("driver.jl")

function run_experiment(params)
  outfile = datadir("BiotBrinkmanVorticityPreconditionerTests",gitdescribe(),
                    savename("BiotBrinkmanVorticityPreconditionerTests",params,"bson"))
  if isfile(outfile)
    println("$outfile (done already)")
    return nothing
  end
  println("$outfile (running)")
  @unpack nk,order,μ,λ,ν,κ,α,c_0,prec_variant = params
  dict = Driver.driver(nk, order, μ, λ, ν, κ, α, c_0, prec_variant)
  merge!(dict,params)
  # @tagsave: add current git commit to dict and then save.
  # "replace_strings_by_symbols" is required to ensure that
  # the dictionary is not of type Dict{Any} but of type
  # Dict{Symbol}
  @tagsave(outfile,replace_strings_by_symbols(dict))
  println(" (done)")
end

# Run all parameter value combinations
dicts=generate_param_dicts()
for params in dicts
  GC.gc()
  run_experiment(params)
end
