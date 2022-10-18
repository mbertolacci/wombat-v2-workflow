#include <algorithm>
#include <map>
#include <tuple>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector to_basis_vector(
    IntegerVector inventory,
    IntegerVector component,
    IntegerVector region,
    IntegerVector month,
    IntegerVector inventory_levels,
    IntegerVector component_levels,
    IntegerVector region_levels,
    IntegerVector month_levels,
    CharacterVector levels
) {
    std::map<
        std::tuple<int, int, int, int>,
        int
    > partsToLevels;

    for (int i = 0; i < levels.size(); ++i) {
        partsToLevels[std::make_tuple(
            inventory_levels[i],
            component_levels[i],
            region_levels[i],
            month_levels[i]
        )] = i + 1;
    }

    IntegerVector output(inventory.size());
    #pragma omp parallel for
    for (int i = 0; i < inventory.size(); ++i) {
        output[i] = partsToLevels[std::make_tuple(
            inventory[i],
            component[i],
            region[i],
            month[i]
        )];
    }
    output.attr("levels") = levels;
    output.attr("class") = "factor";

    return output;
}

// [[Rcpp::export]]
LogicalVector sensitivities_already_done(
    IntegerVector observation_id,
    IntegerVector basis_vector,
    IntegerVector complete_observation_id,
    IntegerVector complete_basis_vector
) {
    LogicalVector output(observation_id.size());
    #pragma omp parallel for
    for (int i = 0; i < observation_id.size(); ++i) {
        auto parts = std::equal_range(
            complete_observation_id.begin(),
            complete_observation_id.end(),
            observation_id[i]
        );
        int start = std::distance(complete_observation_id.begin(), parts.first);
        int end = std::distance(complete_observation_id.begin(), parts.second);
        output[i] = std::binary_search(
            complete_basis_vector.begin() + start,
            complete_basis_vector.begin() + end,
            basis_vector[i]
        );
    }
    return output;
}

// [[Rcpp::export]]
NumericMatrix sparse_triplet_multiply(
    IntegerVector i,
    IntegerVector j,
    NumericVector x,
    LogicalVector zero,
    int n,
    NumericMatrix b
) {
    NumericMatrix output(n, b.ncol());
    std::fill(output.begin(), output.end(), 0);
    #pragma omp parallel for
    for (int k = 0; k < b.ncol(); ++k) {
        for (int index = 0; index < i.size(); ++index) {
            if (zero[index]) continue;
            output(
                i[index] - 1,
                k
            ) += x[index] * b(j[index] - 1, k);
        }
    }
    return output;
}

// [[Rcpp::export]]
IntegerVector update_levels(
    IntegerVector x,
    CharacterVector levels
) {
    CharacterVector oldLevels = x.attr("levels");
    std::map<String,int> levelToIndex;
    for (int i = 0; i < levels.size(); ++i) {
        levelToIndex[levels[i]] = i;
    }
    std::vector<int> oldToNew(oldLevels.size());
    for (int i = 0; i < oldLevels.size(); ++i) {
        oldToNew[i] = levelToIndex[oldLevels[i]];
    }
    IntegerVector output(x.size());
    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i) {
        if (IntegerVector::is_na(x[i])) {
            output[i] = NA_INTEGER;
        } else {
            output[i] = oldToNew[x[i] - 1] + 1;
        }
    }
    output.attr("levels") = levels;
    output.attr("class") = "factor";
    return output;
}

// [[Rcpp::export]]
LogicalVector is_level_in(
    IntegerVector x,
    CharacterVector levels
) {
    CharacterVector xLevels = x.attr("levels");
    std::set<String> levelSet;
    for (int i = 0; i < levels.size(); ++i) {
        levelSet.insert(levels[i]);
    }
    std::vector<bool> keepLevel(xLevels.size());
    for (int i = 0; i < xLevels.size(); ++i) {
        keepLevel[i] = levelSet.count(xLevels[i]) == 1;
    }
    LogicalVector output(x.size());
    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i) {
        output[i] = keepLevel[x[i] - 1];
    }
    return output;
}