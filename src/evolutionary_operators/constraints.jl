function perimeter_constraint(individual::VoronoiIndividual,perimiter_constraint::Float64,exponent = 1)
    return sum(individual.perimeters) < perimiter_constraint ? true : false
end

function area_constraint(individual::VoronoiIndividual,area_constraint::Float64,exponent = 1)
    return sum(individual.areas)^exponent < area_constraint ? true : false
end

function number_constraint(individual::VoronoiIndividual,number_constraint::Float64,exponent = 1)
    return length(individual.genome)  ==  number_constraint ? true : false
end
