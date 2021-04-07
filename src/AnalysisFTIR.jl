module AnalysisFTIR

using MultiLayerNFRHT, Plots; pyplot()
using CSV, DataFrames, LsqFit

export load_data,convertto_freq,convertto_wavel,
       unnormalize,mynormalize,plot_data,reflectivity_data,
       emissivity_data,comparison,emissivity_model,
       save_data, fit_data,
       window_average, savitsky_golay


const dfltoptions = Dict("Au layer thickness (nm)" => 1,
                         "Electron mean free path (nm)" => 2 ,
                         "Si layer thickness (nm)" => 3 )


" Fitting options "
function option(opt :: Symbol)
    return options[opt]
end


" Load the data and put into an array "
function load_data(filename)
    df    = DataFrame(CSV.File(filename))
    mdata = convert(Array,df)
    return mdata
end


" Convert cm^{-1} to frequency in rad/s "
convertto_freq(wavenum) = 1e2*2.0*pi*c0*wavenum

" Convert cm^{-1} to wavelength in microns "
convertto_wavel(wavenum) = 1e4/wavenum

" Remove the reference from the experimental data"
unnormalize(data,refdata) = data[:,2].*refdata[:,2]

" Normalize the experimental data"
mynormalize(data,normdata) = data[:,2]./normdata[:,2]

" Plot reflectivity/emissivity as a function of frequency"
function plot_data(data1,data2)
    plot(data1[:,1],data1[:,2]; xscale= :log10,ylim=(0,1))
    plot!(data2[:,1],data2[:,2]; xscale= :log10,ylim=(0,1))
end

"""
Extract reflectivity with a given reference where srcdata, initrefdata and refdata are filenames.
initrefdata is the data that was defined as reference data during the measurements.
All measurements are normalized by this dataset. All source data have to be multiplied by this set.
"""
function reflectivity_data(srcdata,initrefdata,refdata)
    data1      = load_data(srcdata)
    initrefdat = load_data(initrefdata)
    refdat     = load_data(refdata)

    data2      = zeros(refdat)
    data2[:,1] = convertto_freq.(data1[:,1])
    data2[:,2] = mynormalize(data1,refdat) # since both datasets are devided by the same signal, unnormalize is not necessary
    return data2
end

" Get the emissivity from the reflectivity measurement "
function emissivity_data(srcdata,initrefdata,refdata)
    data      = reflectivity_data(srcdata,initrefdata,refdata)
    data[:,2].= 1.0 .- data[:,2]
    return data
end

" Compare simulation and measurement "
function comparison(srcdata,initrefdata,refdata,substrate :: MultiLayer; save = false, outputfilename = "")
    emdata   = emissivity_data(srcdata,initrefdata,refdata)
    emw      = copy(emdata)
    emw[:,2] = emissivity_kx_w.([substrate],0.0, emdata[:,1]).*0.5
    if save == true
        save_data(emdata,emw,outputfilename)
    end
    plot_data(emdata,emw)
end


" Fit model to experimental data "
function fit_data(emissivity_model, xdata, ydata, p; save = false, outputfilename = "" ,opt=dfltoptions)
    fit     = curve_fit(emissivity_model,xdata,ydata,p)
    fit_dat = zeros(xdata)
    fit_dat = emissivity_model(xdata, fit.param)
    df_fit = value_of_fit(fit,opt)
    if save == true
        save_data([xdata  ydata] ,[xdata fit_dat], outputfilename)
        writetable(outputfilename*"_parameters",df_fit)
    end
    data = [xdata  ydata  fit_dat]
    println(df_fit)
    plot_data([xdata ydata],[xdata fit_dat])
end

function value_of_fit(fit, dict)
    df = DataFrame()
    df[:Parameters] = [collect(keys(dict)) ; "Sum residuals"]
    df[:Values]     = [[fit.param[i] for i in collect(values(dict))] ; sum(fit.resid.^2)]
    return df
end


" Save experimental and simulated emissivity to a *.dat file"
function save_data(data_exp,data_sim,outputfilename)
    df        = DataFrame()
    df[:freq] = data_exp[:,1]
    df[:Data] = data_exp[:,2]
    df[:Fit]  = data_sim[:,2]
    writetable(outputfilename, df,separator = ' ')
end

" Sliding window average "
function window_average(data :: AbstractArray, window :: Integer)
    data_av = zeros(length(data))
    for i=1:length(data)
        data_av[i] = mean(data[i:min(i+window-1,end)])
    end
    return data_av
end

" Savitsky-Golay filter for smoothing data "
function savitsky_golay(x::Vector, windowSize::Integer, polyOrder::Integer; deriv::Integer=0)

    #Some error checking
    @assert isodd(windowSize) "Window size must be an odd integer."
    @assert polyOrder < windowSize "Polynomial order must me less than window size."

    halfWindow = Int((windowSize-1)/2)

    #Setup the S matrix of basis vectors.
    S = zeros(windowSize, polyOrder+1)
    for ct = 0:polyOrder
	    S[:,ct+1] = collect(-halfWindow:halfWindow).^(ct)
    end

    #Compute the filter coefficients for all orders
    #From the scipy code it seems pinv(S) and taking rows should be enough
    G = S*pinv(S'*S)

    #Slice out the derivative order we want
    filterCoeffs = G[:,deriv+1] * factorial(deriv)

    #Pad the signal with the endpoints and convolve with filter
    paddedX = [x[1]*ones(halfWindow); x ; x[end]*ones(halfWindow)]
    y = conv(filterCoeffs[end:-1:1], paddedX)

    #Return the valid midsection
    return y[2*halfWindow+1:end-2*halfWindow]

end


end # module
