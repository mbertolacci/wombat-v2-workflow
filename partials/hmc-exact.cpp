#include <RcppEigen.h>
#include <map>

// [[Rcpp::depends(RcppEigen)]]

using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::NumericVector;
using Rcpp::S4;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
typedef Eigen::SparseMatrix<double> MatrixSd;
typedef Eigen::Map<MatrixSd> MapSd;
typedef Eigen::Map<MatrixXd> MapXd;

// [[Rcpp::export]]
NumericVector sampleHmcConstrained(
    VectorXd x0,
    VectorXd mu,
    MatrixXd R,
    MatrixSd F,
    VectorXd g,
    double totalTime,
    double tolerance = 1e-12,
    int bounceLimit = 1000,
    int nSamples = 1,
    bool debug = false
) {
    VectorXd gPrime = g + F * mu;
    x0 = R.triangularView<Eigen::Upper>() * (x0 - mu);

    VectorXd c = F * R.triangularView<Eigen::Upper>().solve(x0) + gPrime;
    if ((c.array() < 0.0).any()) {
        Rcpp::stop("x0 does not satisfy constraints");
    }

    std::map<int,VectorXd> fPrimeRowCache;
    std::map<int,double> fPrimeNormCache;

    VectorXd x = x0;
    MatrixXd samples(nSamples, x0.size());

    Rcpp::Rcout.precision(4);

    for (int iteration = 0; iteration < nSamples; ++iteration) {
        int lastBounceIndex = -1;
        double timeTaken = 0;
        VectorXd v = Rcpp::as<VectorXd>(Rcpp::rnorm(x0.size()));

        int i;
        for (i = 0; i < bounceLimit; ++i) {
            if (i % 10 == 0) Rcpp::checkUserInterrupt();
            if (debug) {
                Rcpp::Rcout
                    << "\rsample = " << (iteration + 1)
                    << ", bounce = " << i
                    << ", timeTaken = " << std::fixed << timeTaken
                    << "/" << std::fixed << totalTime
                    << ", cache size = " << fPrimeRowCache.size();
            }

            VectorXd a = v;
            VectorXd b = x;
            ArrayXd f1 = (F * R.triangularView<Eigen::Upper>().solve(a)).array();
            ArrayXd f2 = (F * R.triangularView<Eigen::Upper>().solve(b)).array();

            ArrayXd U = (f1.square() + f2.square()).sqrt();
            ArrayXd phi = (-f1).binaryExpr(f2, [] (double a, double b) { return std::atan2(a, b); });

            auto canHit = (gPrime.array() / U).abs() <= 1;
            double deltaTime = totalTime;
            if (canHit.any()) {
                ArrayXd hitTime = -phi + (-gPrime.array() / U).acos();

                if (lastBounceIndex > -1 && canHit[lastBounceIndex]) {
                    if (
                        std::abs(hitTime[lastBounceIndex]) < tolerance
                        || std::abs(hitTime[lastBounceIndex] - 2 * M_PI) < tolerance
                    ) {
                      hitTime[lastBounceIndex] = std::numeric_limits<double>::max();
                    }
                }

                double leastTime = std::numeric_limits<double>::max();
                int leastTimeIndex = -1;
                for (int j = 0; j < canHit.size(); ++j) {
                    if (canHit[j] && hitTime[j] < leastTime) {
                        leastTimeIndex = j;
                        leastTime = hitTime[j];
                    }
                }
                deltaTime = leastTime;
                lastBounceIndex = leastTimeIndex;
            }

            timeTaken += deltaTime;
            if (timeTaken >= totalTime) {
                deltaTime -= timeTaken - totalTime;
            }

            x = a * sin(deltaTime) + b * cos(deltaTime);
            v = a * cos(deltaTime) - b * sin(deltaTime);

            if (timeTaken >= totalTime) {
                break;
            }

            // compute reflected velocity
            if (fPrimeRowCache.count(lastBounceIndex) == 0) {
                VectorXd fRowCurrent = F.row(lastBounceIndex);
                fPrimeRowCache[lastBounceIndex] = R.triangularView<Eigen::Upper>()
                    .transpose()
                    .solve(fRowCurrent);
                fPrimeNormCache[lastBounceIndex] = fPrimeRowCache[lastBounceIndex]
                    .squaredNorm();
            }
            double qReflection = (
                fPrimeRowCache[lastBounceIndex].array() * v.array()
            ).sum() / fPrimeNormCache[lastBounceIndex];
            v = v - 2 * qReflection * fPrimeRowCache[lastBounceIndex];
        }
        if (debug) {
            Rcpp::Rcout << "\n";
        }

        if (i == bounceLimit) {
            Rcpp::stop("Bounce limit exceeded");
        }
        samples.row(iteration) = (
            mu + R.triangularView<Eigen::Upper>().solve(x)
        ).transpose();
    }

    if (debug) {
        Rcpp::Rcout << "final cache size = " << fPrimeRowCache.size() << "\n";
    }

    NumericVector output;
    if (nSamples == 1) {
        output = Rcpp::wrap(samples.row(0).transpose());
    } else {
        output = Rcpp::wrap(samples);
    }

    return output;
}
