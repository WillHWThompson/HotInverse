function perimeter_constraint(individual::VoronoiIndividual,perimiter_constraint::Float64,exponent::Float64 = 1)
    ind_sum = sum(individual.perimeters)
    bool_val =  ind_sum < perimiter_constraint ? true : false
    @debug "constraint_sum $ind_sum, constraint_value: $perimiter_constraint" #$(sum(individual.areas)"
    return bool_val
end

function area_constraint(individual::VoronoiIndividual,area_constraint::Float64,exponent::Float64 = 1)
    ind_sum = sum(individual.areas)^exponent 
    bool_val = ind_sum < area_constraint ? true : false
    @debug "constraint_sum $ind_sum, constraint_value: $area_constraint" #$(sum(individual.areas)"
    return bool_val
end

function number_constraint(individual::VoronoiIndividual,number_constraint::Float64,exponent::Float64 = 1)
    ind_sum = length(individual.genome)
    bool_val = ind_sum  ==  number_constraint ? true : false
    @debug "constraint_sum $ind_sum, constraint_value: $number_constraint" #$(sum(individual.areas)"
    return bool_val
end
