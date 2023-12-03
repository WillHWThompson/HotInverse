function perimeter_constraint(individual::VoronoiIndividual,perimiter_constraint::Float64,exponent = 1;return_val = false)
    bool_val = sum(individual.perimeters) < perimiter_constraint ? true : false
    if return_val
        return bool_val,sum(individual.perimeters),perimiter_constraint 
    else
        return bool_val
    end
end

function area_constraint(individual::VoronoiIndividual,area_constraint::Float64,exponent = 1;return_val = false)
    bool_val = sum(individual.areas)^exponent < area_constraint ? true : false
    if return_val
        return bool_val,sum(individual.areas)#,area_constraint 
    else
        return bool_val
    end
end

function number_constraint(individual::VoronoiIndividual,number_constraint::Float64,exponent = 1;return_val = false)
    bool_val = length(individual.genome)  ==  number_constraint ? true : false
    if return_val
        return bool_val,length(individual.genome),number_constraint 
    else
        return bool_val
    end
end
