// Reference: Hoots, Felix R. and Roehrich, Ronald L.,
// "SPACETRACK REPORT NO. 3: Models for Propagation of NORAD Element Sets,"
// December 1980.

#ifndef SGP4_H
#define SGP4_H

#include <string>

namespace sgp4
{
    class Satellite
    {
    public:

        // Create satellite from two-line element set.
        Satellite(
            const std::string& line1 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    87",
            const std::string& line2 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  1058");

        // Return 5-digit catalog number.
        long get_catalog_number() const;

        // Return epoch as Julian day number (UTC).
        double get_epoch() const;

        // Propagate to given time (seconds since epoch) and set ECI state to
        // {x, y, z, xdot, ydot, zdot} (true equator mean equinox of epoch,
        // meters and meters/second).
        void coast(double tsince, double *eci) const;

        // Convert state from ECI to EFG at given Julian day number (UT1).
        void eci_to_efg(double jd, const double *eci, double *efg) const;

    protected:

        // Element set fields
        long catalog_number;
        char classification;
        std::string designator;
        double epoch;
        double xndt2o;
        double xndd6o;
        double bstar;
        char ephemeris_type;
        long elset_number;
        double xincl;
        double xnodeo;
        double eo;
        double omegao;
        double xmo;
        double xno;
        long rev_number;

        bool isimp;
        double xmdot, omgdot, xnodot, xnodcf, c1, c4, t2cof, omgcof, xmcof,
            eta, delmo, d2, d3, d4, c5, sinmo, t3cof, t4cof, t5cof, aodp,
            xnodp, xlcof, aycof, x3thm1, x1mth2, x7thm1, cosio, sinio;
    };
} // namespace sgp4

#endif // SGP4_H
