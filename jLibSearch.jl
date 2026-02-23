# Enable using packages
using Pkg

# Install and update packages
Pkg.update() # Ensure everything is up-to-date
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Glob")
Pkg.add("XLSX")
Pkg.add("Statistics")
Pkg.add("LinearAlgebra")
Pkg.add("PlotlyJS")

# Import packages
using CSV
using DataFrames
using XLSX
using Statistics
using LinearAlgebra
using Dates
using PlotlyJS



#######################################################################
############ Functions for library search & visualization #############
#######################################################################



# A function to merge RI values to the leftmost non-missing value 
# and pair them with corresponding RI tolerance values
function mergeToLeft(ri, RIwin)
    RIval = fill(NaN, size(ri, 1))
    RItol = fill(NaN, size(ri, 1))
    RI_Type = fill("", size(ri, 1))
    for r = 1:size(ri, 1)
        if ismissing(ri[r, 1])
            if !ismissing(ri[r, 2])
                RIval[r] = ri[r, 2]
                RItol[r] = RIwin[2]
                RI_Type[r] = "Std_NP"
            elseif !ismissing(ri[r, 3])
                RIval[r] = ri[r, 3]
                RItol[r] = RIwin[3]
                RI_Type[r] = "AI_Semi-Std_NP"
            end
        else
            RIval[r] = ri[r, 1]
            RItol[r] = RIwin[1]
            RI_Type[r] = "Semi-Std_NP"
        end
    end
    RIval = float.(RIval)
    RIval[ismissing.(RIval)] .= NaN

    return RIval, RItol, RI_Type
end



# A function to split spectrum and turn it into a 2D array
function getSpec(sp)
    # Remove trailing space if present
    if sp[end] == ' '
        sp = sp[1:end-1]
    end

    # Check if spectrum is Sync 2D type
    if occursin(r"\(\d+\|\d+\.\d+\)", sp)
        # Convert Sync 2D type spectra to LECO ChromaTOF type spectra
        sp = replace(sp, r"\((\d+)\|(\d+\.\d+)\)" => s"\1:\2 ")
        sp = replace(sp, r"\s*\)\s*\(" => " ")
        sp = replace(sp, r"^\(|\)$" => "")
    end

    # Remove trailing space if present
    if sp[end] == ' '
        sp = sp[1:end-1]
    end

    # Split the string by commas & spaces (Guineu) or spaces (LECO ChromaTOF ) and then by colons
    # For Guineu spectra
    if sp[1] == '['
        sp = sp[2:end-1]
        sp = split(sp, " , ")
    # For LECO ChromaTOF spectra
    else
        sp = split(sp, " ")
    end
    sp = split.(sp, ":")
    # Parse the split strings into Float64 and reshape into a 2D array
    sp = parse.(Float64, stack(sp, dims = 1))

    return sp
end



# A function to convert Sync 2D type spectrum into LECO ChromaTOF type spectrum
function SpecSync2DtoLECO(sp)
    # Remove trailing space if present
    if sp[end] == ' '
        sp = sp[1:end-1]
    end

    # Check if spectrum is Sync 2D type
    if occursin(r"\(\d+\|\d+\.\d+\)", sp)
        # Convert Sync 2D type spectra to LECO ChromaTOF type spectra
        sp = replace(sp, r"\((\d+)\|(\d+\.\d+)\)" => s"\1:\2 ")
        sp = replace(sp, r"\s*\)\s*\(" => " ")
        sp = replace(sp, r"^\(|\)$" => "")
    end

    # Remove trailing space if present
    if sp[end] == ' '
        sp = sp[1:end-1]
    end

    # Split the spectrum into pairs
    pairs = split(sp, ' ')
    
    # Find the maximum intensity value for rescaling
    intensities = [parse(Float64, split(pair, ':')[2]) for pair in pairs]
    max_intensity = maximum(intensities)

    # Round values and rescale intensity values
    rounded_rescaled_pairs = map(pair -> begin
        mz, intensity = split(pair, ':')
        mz = round(parse(Int, mz))
        intensity = round(parse(Float64, intensity))
        rescaled_intensity = round(intensity * 9999 / max_intensity)
        "$mz:$rescaled_intensity"
    end, pairs)

    # Join the pairs back into a single string
    sp = join(rounded_rescaled_pairs, ' ')

    return sp
end



# A function to generate a master list of m/z values and their corresponding intensities from two sets of input data
function master_mass_gen(mz_ref, int_ref, mz_user, int_user, mass_tol)
    mz_values = vcat(mz_ref[:], mz_user[:])
    master_vect = zeros(length(vcat(mz_ref[:], mz_user[:])), 3)

    for i = 1:size(master_vect, 1)
        tv1 = abs.(mz_values[i] .- mz_values)
        tv2 = findall(x -> x <= mass_tol/2, tv1)
        master_vect[i, 1] = mean(mz_values[tv2])
        mz_values[tv2] .= 0
    end

    for i = 1:size(master_vect, 1)
        if master_vect[i] > 0
            tv1 = abs.(master_vect[i] .- mz_ref)
            tv2 = findall(x -> x <= mass_tol/2, tv1)
            if length(tv2) > 0
                master_vect[i, 2] = maximum(int_ref[tv2])
                int_ref[tv2] .= 0
            end

            ttv1 = abs.(master_vect[i] .- mz_user)
            ttv2 = findall(x -> x <= mass_tol/2, ttv1)
            if  length(ttv2) > 0
                if maximum(ttv2) <= length(int_user)
                    master_vect[i, 3] = maximum(int_user[ttv2])
                    int_user[ttv2] .= 0
                elseif maximum(ttv2) > length(int_user)
                    tv3 = zeros(size(int_user))
                    for j = 1:length(ttv2)
                        if ttv2[j] > length(int_user)
                            ttv2[j] = 0
                        end

                    end
                    master_vect[i, 3] = maximum(int_user[ttv2[ttv2 .> 0]])
                    int_user[ttv2[ttv2 .> 0]] .= 0
                end
            end
        end
    end

    master_vect_f = master_vect[master_vect[:, 1] .> 0, :]

    return master_vect_f
end



# A function to search library/database to identify unknowns and verify previous identifications
function librarySearch(pathDB, pathFile, RIwin, mz_tol, numfrags = 15, similarity_method::String="DISCO")
    println("*"^150)
    println("Searching library/database... It can take a while... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    time_b4_search = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("-"^150)

    # Load database and file
    xf = XLSX.readxlsx(pathDB)
    sn = XLSX.sheetnames(xf)
    db = DataFrame(XLSX.readtable(pathDB, sn[1]))

    # Create decided RI columns
    RIval, RItol, RI_Type = mergeToLeft(Matrix(db[:,3:12])[:,1:3:end], RIwin)
    db[!, "RImatching"] = RIval
    db[!, "RImatchingTol"] = RItol
    db[!, "RImatchingType"] = RI_Type
    spec = []

    # Iterate candidate spectra    
    if (split(pathFile,".")[end] .== "csv" ) || (split(pathFile,".")[end] .== "txt" )
        df = CSV.read(pathFile, DataFrame)
        if "R.I. calc" in names(df)
            colRI = "R.I. calc"
            colNr = "ID"
        else
            colRI = "Retention Index"
            colNr = "Nr"
        end
    elseif split(pathFile,".")[end] .== "xlsx"
        xf = XLSX.readxlsx(pathFile)
        sn = XLSX.sheetnames(xf)
        df = DataFrame(XLSX.readtable(pathFile, sn[1]))

        # Check for Sync 2D vs. Guineu input file
        if "R.I. calc" in names(df)
            colRI = "R.I. calc"
            colNr = "ID"
        else
            colRI = "RTI"
            colNr = "ID"
        end
    else
        println("File format not reconized.")
        return
    end
    df.Nr = 1:size(df, 1)

    if any(names(df) .== "Spectrum")
        spec = df[:, "Spectrum"]
    else
        println("Spectrum column not found. Looking for: \"Spectrum\".")
        return
    end

    # Convert "Feature #" in "Name" column to "Feature_RT1_RT2"
    for i in 1:size(df, 1)
        if occursin(r"Feature \d+", df[i, "Name"])
            RT1 = round(Int, df[i, "Med RT1 (sec)"])
            RT2 = round(df[i, "Med RT2 (sec)"], digits=3)
            df[i, "Name"] = "Feature_$(RT1)_$(RT2)"
        end
    end
    
    # Setup empty dataframe
    identAll = DataFrame(ID = 0, DB_ID = 0, MinRI = 0., MaxRI = 0., RI = 0., DeltaRI = 0., DB_RI = 0., DB_RI_Type = "", Dot_Score = 0., CL = "", Name = "", DB_Name = "", DB_Formula = "", DB_InChIKey = "", DB_CAS = "", Spectrum = "", DB_Spectrum = "", )

    for i = 1:size(df, 1)
        println(i)
        sp = getSpec(spec[i])

        # Select candidates
        minRI = fill(df[i, colRI], size(db, 1)) - db[!, "RImatchingTol"]
        maxRI = fill(df[i, colRI],size(db, 1)) + db[!, "RImatchingTol"]

        ind_s = findall(minRI .<= db[!, "RImatching"] .<= maxRI)
        if isempty(ind_s)
            continue
        end
        ident = DataFrame(ID = 0, DB_ID = 0, MinRI = 0., MaxRI = 0., RI = 0., DeltaRI = 0., DB_RI = 0., DB_RI_Type = "", Dot_Score = 0., CL = "", Name = "", DB_Name = "", DB_Formula = "", DB_InChIKey = "", DB_CAS = "", Spectrum = "", DB_Spectrum = "", )

        for s in ind_s
            # Calculate dot product
            spdb = getSpec(db[s, "Spectrum"])
            master_vect_f = master_mass_gen(sp[:, 1], sp[:, 2], spdb[:, 1], spdb[:, 2], mz_tol)
            
            if size(master_vect_f, 1) .> numfrags
                master_vect_fu = master_vect_f[sortperm(vec(sum(master_vect_f[:, 3],dims = 2)), rev = true), :][1:numfrags, :]
                master_vect_fr = master_vect_f[sortperm(vec(sum(master_vect_f[:, 2],dims = 2)), rev = true), :][1:numfrags, :]
                master_vect_f = [master_vect_fu; master_vect_fr]
                master_vect_f = unique(master_vect_f, dims = 1)
            end
            # Mean center each vector or not based on method
            if similarity_method == "DISCO"
                # Mean center intensities
                master_vect_f[:, 2:3] = master_vect_f[:, 2:3] .- mean(master_vect_f[:, 2:3], dims = 1)
                p = abs(dot(master_vect_f[:, 2], master_vect_f[:, 3])/(norm(master_vect_f[:, 2]*norm(master_vect_f[:, 3]))))
            elseif similarity_method == "NDP"
                # Do not mean-center; just normalize directly
                p = dot(master_vect_f[:, 2], master_vect_f[:, 3]) / (norm(master_vect_f[:, 2]) * norm(master_vect_f[:, 3]))
            else
                error("Unknown similarity method: $similarity_method. Use \"DISCO\" or \"NDP\".")
            end
            p = round(p, digits=2) # Round "Dot_Score" to 2 decimal places

            # Safe column value getter to account for columns that do not exist
            getval(row, col, default="") = haskey(row, col) ? row[col] : default

            ident = push!(ident, (
                getval(df[i, :], colNr, 0),
                getval(db[s, :], "#", 0),
                minRI[s],
                maxRI[s],
                getval(df[i, :], colRI, 0.0),
                getval(df[i, :], colRI, 0.0) - getval(db[s, :], "RImatching", 0.0),
                getval(db[s, :], "RImatching", 0.0),
                getval(db[s, :], "RImatchingType", ""),
                p,
                getval(df[i, :], "CL Sync2D", ""),
                getval(df[i, :], "Name", ""),
                getval(db[s, :], "Name", ""),
                getval(db[s, :], "Formula", ""),
                getval(db[s, :], "InChIKey", ""),
                string(getval(db[s, :], "CAS Number", "")),
                SpecSync2DtoLECO(spec[i]),
                SpecSync2DtoLECO(getval(db[s, :], "Spectrum", ""))),
                promote = true
            )
        end

        # Only keep top 3 matches
        # Define the priority for "DB_RI_Type" values
        priority = Dict("Semi-Std_NP" => 1, "Std_NP" => 2, "AI_Semi-Std_NP" => 3, "" => 4, missing => 4)
        # Sort ident based on "DB_RI_Type" priority and "Dot_Score" values
        ident = ident[sortperm(ident[:, "DB_RI_Type"], by = x -> priority[x]), :]
        ident = ident[sortperm(ident[:, "Dot_Score"], rev = true), :]
        if size(ident, 1) > 3
            ident = ident[1:3, :]
        end
        # Add to all identifications
        identAll = [identAll; ident]
    end

    identAll = identAll[2:end, :]

    if (split(pathFile,".")[end] .== "csv" ) || (split(pathFile,".")[end] .== "txt" )
        CSV.write(pathFile[1:end-4]*"_LibSearch.csv", identAll)
    elseif split(pathFile,".")[end] .== "xlsx"
        CSV.write(pathFile[1:end-5]*"_LibSearch.csv", identAll)
    end
    println("-"^150)
    println("Done library searching.")
    println("Start time: ", time_b4_search)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)    
end



# Helper functions to cleanup compound names
# Keep only printable ASCII (32..126), drop ALL control chars like \r \n \t
printable_ascii(c) = (' ' <= c <= '~')
function ascii_printable_only(s::AbstractString)
    io = IOBuffer()
    for c in String(s)
        if printable_ascii(c)
            write(io, c)
        end
    end
    return String(take!(io))
end

# Targeted replacements for common unicode science symbols BEFORE ASCII-only filter
function replace_common_unicode(s::AbstractString)
    str = String(s)
    # Common mappings
    pairs = [
        ("β","beta"), ("α","alpha"), ("γ","gamma"), ("δ","delta"),
        ("µ","mu"),   ("ß","beta"), ("a-","alpha-"), ("-a-","-alpha-"),
        ("–","-"), ("—","-"),
        ("’","'"), ("‘","'"), ("“","\""), ("”","\"")
    ]
    for (sym, repl) in pairs
        str = replace(str, string(sym) => repl)
    end
    return str
end

# Final sanitizer: map known unicodes, strip control chars, drop non-ASCII
sanitize_for_json(s::AbstractString) = ascii_printable_only(replace_common_unicode(s))

# A helper function to wrap long text every n characters (works safely after ASCII sanitizing)
function wrap_text(s::AbstractString, n::Int=20)
    s2 = sanitize_for_json(String(s))   # explicitly convert to standard String
    return isempty(s2) ? s2 : join([s2[i:min(i+n-1,end)] for i in 1:n:length(s2)], "<br>")
end



# A function to create head-to-tail graphs for user vs. library spectra
function libraryVisualization(pathFile, index)
    time_b4_vis = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Creating graphs... It can take a while... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
   
    df = CSV.read(pathFile, DataFrame)
    df[!, "Graph"] = fill("NA", size(df, 1))
    # Debugging: delete rows where 'Spec' column contains missing values
    dropmissing!(df, :Spectrum)

    if typeof(index) == String
        index = 1:size(df, 1)
    end

    # Create a folder for the graphs
    pathout = pathFile[1:end-4]
    if !isdir(pathout)
        mkdir(pathout)
    end

    for i in index
        us = getSpec(df[i, "Spectrum"])
        db = getSpec(df[i, "DB_Spectrum"])

        # Normalize intensities; guard against zero/NaN maxima (avoid NaNs in JSON)
        usmax = maximum(us[:,2])
        dbmax = maximum(db[:,2])
        us_y = (isfinite(usmax) && usmax > 0) ? (us[:,2] ./ usmax) .* 100 : zeros(size(us,1))
        db_y = (isfinite(dbmax) && dbmax > 0) ? (db[:,2] ./ dbmax) .* -100 : zeros(size(db,1))

        # Base head-to-tail spectrum
        trace_user = bar(
            x=us[:,1], y=us_y, name="User spectrum",
            marker_color="rgb(54,96,146)", width=1.0, offset=0.0
        )

        trace_db = bar(
            x=db[:,1], y=db_y, name="Database spectrum",
            marker_color="rgb(150,54,52)", width=1.0, offset=0.0
        )

        # Annotate peaks: user spectrum
        ann = []
        bit = ones(size(us,1))
        while any(bit .== 1)
            ind = argmax(us[:,2] .* bit)
            push!(ann, attr(x=us[ind,1], y=us_y[ind]+8,
                            text=string(Int(us[ind,1])), showarrow=false, font=attr(size=9)))
            bit[us[ind,1]-5 .< us[:,1] .<= us[ind,1]+5] .= 0
        end

        # Annotate peaks: database spectrum
        bit = ones(size(db,1))
        while any(bit .== 1)
            ind = argmax(db[:,2] .* bit)
            push!(ann, attr(x=db[ind,1], y=db_y[ind]-8,
                            text=string(Int(db[ind,1])), showarrow=false, font=attr(size=9)))
            bit[db[ind,1]-5 .< db[:,1] .<= db[ind,1]+5] .= 0
        end

        # Sanitize feature names
        user_name = wrap_text(string(df[i,"Name"]), 20)
        db_name   = wrap_text(string(df[i,"DB_Name"]),   20)

        # Feature name annotations outside the plot area
        push!(ann, attr(
            x=1.005, # Slightly right of plot area (in paper coordinates)
            y=0.75, # Top spectrum, relative to 0-1 paper coordinates
            xref="paper",
            yref="paper",
            text="User: "*user_name,
            showarrow=false,
            xanchor="left",
            yanchor="middle",
            textangle=-90,
            font=attr(size=14, color="rgb(54,96,146)")
        ))
        push!(ann, attr(
            x=1.005, # Slightly right of plot area
            y=0.25, # Bottom spectrum, relative to 0-1 paper coordinates
            xref="paper",
            yref="paper",
            text="DB: "*db_name,
            showarrow=false,
            xanchor="left",
            yanchor="middle",
            textangle=-90,
            font=attr(size=14, color="rgb(150,54,52)")
        ))

        # Layout with secondary y-axis for feature names
        layout = Layout(
            title=attr(text="Head-to-tail spectrum comparison", x=0.5), # Centered title
            xaxis=attr(title="m/z"),
            yaxis=attr(title="Relative intensity", range=[-110,110], zeroline=true, zerolinecolor="black"),
            bargap=0.05,
            showlegend=false,
            annotations=ann,
            margin=attr(r=100)  # Increase r for more white space on the right
        )

        fig = Plot([trace_user, trace_db], layout) 
        
        # Save as PNG
        fname = "$(i)_User_$(df[i,"ID"])=DB_$(df[i,"DB_ID"]).png"
        fpath = joinpath(pathout, fname)
        try
            savefig(fig, fpath)
            df[i,"Graph"] = "file:///" * fpath
        catch err
            # Retry with annotations removed (if text still problematic)
            try
                layout_no_ann = deepcopy(layout); layout_no_ann[:annotations] = Any[]
                fig2 = Plot([trace_user, trace_db], layout_no_ann)
                savefig(fig2, fpath)
                df[i,"Graph"] = "file:///" * fpath
                @warn "Saved without name annotations for entry $(i) due to text encoding."
            catch err2
                # Final fallback: save HTML (always works) so you at least keep a visual
                htmlpath = replace(fpath, ".png" => ".html")
                PlotlyJS.savehtml(fig, htmlpath)
                df[i,"Graph"] = "file:///" * htmlpath
                @warn "PNG export failed for entry $(i). Saved HTML instead at $(htmlpath)."
            end
        end
    end
    
    CSV.write(pathFile, df)
    println("\n", "Done creating graphs.", "\n")
    println("Start time: ", time_b4_vis)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
    println("*"^150)
end



#######################################################################
############ Execution of library search & visualization ##############
#######################################################################



# Library search
# Specify library and aligned files
println("Specifying library/database file... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
pathDB = "F:\\Projects\\NIST23_Final_Excel_NoDups_from_mz45_Rescaled_NoPolarRI.xlsx"    # pathDB = "F:\\UMU\\INQUIRE\\Acquired_Data\\PEG_Data\\NIST23_Excel_Short_NoDups.xlsx"
println("Done specifying library/database. ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")

println("Specifying aligned file... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
pathFile = "F:\\Projects\\INQUIRE\\NIST_Dust_Trial_GCxGC\\230727_Three-City_Dust_Trial\\Dust_Trial_Paper\\Di-n-octyl_phthalate_Julia_LibSearch\\251009_Di-n-octyl_phthalate_Target_Search.xlsx"
println("Done specifying aligned file. ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")

# Run library search
RIwin = [50, 60, 100]   # RI tolerance windows for semi-standard non-polar RI, standard non-polar RI, and AI generated RI
mz_tol = 0.1    # m/z tolerance for matching spectra
numfrags = 50       # Minimum number of highest intensity fragments used from each spectrum (default = 15); can be left out for running the function
similarity_method = "DISCO"    # Spectral similarity method: "DISCO" (DIstance & Spectrum Correlation Optimization) or "NDP" (Normalized Dot Product)
librarySearch(pathDB, pathFile, RIwin, mz_tol, numfrags, similarity_method)



# Library search visualization
println("Specifying file for visualization... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
pathFile = "F:\\Projects\\INQUIRE\\NIST_Dust_Trial_GCxGC\\230727_Three-City_Dust_Trial\\Dust_Trial_Paper\\Di-n-octyl_phthalate_Julia_LibSearch\\251009_Di-n-octyl_phthalate_Target_Search_LibSearch.csv"
println("Done specifying file for visualization... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")
index = "all"       # Which features should be plotted: "all" plots every entry, a vector of numbers [1,2,6] plots only those indices, and collect(1:100) plots only indices/entries from 1 to 100
libraryVisualization(pathFile, index)