#include "sgp4.h"
#include <cmath>
#include <cstdlib>

namespace // unnamed namespace
{
    const double E6A = 1.e-6;
    const double QO = 120.0;
    const double SO = 78.0;
    const double TOTHRD = 2. / 3.;
    const double TWOPI = 6.283185307179586;
    const double DE2RA = TWOPI / 360;
    const double XJ2 = 1.082616e-3;
    const double XJ3 = -.253881e-5;
    const double XJ4 = -1.65597e-6;
    const double XKE = .743669161e-1;
    const double XKMPER = 6378.135;
    const double XMNPDA = 1440.;
    const double AE = 1.;

    const double CK2 = .5*XJ2*std::pow(AE, 2.);
    const double CK4 = -.375*XJ4*std::pow(AE, 4.);
    const double QOMS2T = std::pow((QO - SO)*AE / XKMPER, 4.);
    const double S = AE*(1. + SO / XKMPER);

    const double POS_UNITS = XKMPER / AE * 1e3;
    const double VEL_UNITS = POS_UNITS * XMNPDA / 86400;
    const double OMEGA = 0.72921151467e-4;
} // unnamed namespace

sgp4::Satellite::Satellite(const std::string& line1, const std::string& line2)
{
    // Parse line 1 of element set.
    catalog_number = std::strtol(line1.substr(3, 4).c_str(), 0, 10);
    const char *alpha5 = "0123456789ABCDEFGHJKLMNPQRSTUVWXYZ";
    for (int c = 0; c < 34; ++c)
    {
        if (line1[2] == alpha5[c])
        {
            catalog_number += c * 10000;
            break;
        }
    }
    classification = line1[7];
    designator = line1.substr(9, 8);
    long year = 1900 + std::strtol(line1.substr(18, 2).c_str(), 0, 10);
    if (year < 1957)
    {
        year += 100;
    }
    year += 4799;
    epoch = 365 * year + year / 4 - year / 100 + year / 400 - 31739.5 +
        std::strtod(line1.substr(20, 12).c_str(), 0);
    xndt2o = std::strtod(line1.substr(33, 10).c_str(), 0)*
        (TWOPI / XMNPDA / XMNPDA);
    xndd6o = std::strtod(line1.substr(44, 6).c_str(), 0) * std::pow(10., -5 +
        std::strtod(line1.substr(50, 2).c_str(), 0))*
        (TWOPI / XMNPDA / XMNPDA / XMNPDA);
    bstar = std::strtod(line1.substr(53, 6).c_str(), 0) * std::pow(10., -5 +
        std::strtod(line1.substr(59, 2).c_str(), 0)) / AE;
    ephemeris_type = line1[62];
    elset_number = std::strtol(line1.substr(64, 4).c_str(), 0, 10);

    // Parse line 2 of element set.
    xincl = std::strtod(line2.substr(8, 8).c_str(), 0)*DE2RA;
    xnodeo = std::strtod(line2.substr(17, 8).c_str(), 0)*DE2RA;
    eo = std::strtod(line2.substr(26, 7).c_str(), 0)*1e-7;
    omegao = std::strtod(line2.substr(34, 8).c_str(), 0)*DE2RA;
    xmo = std::strtod(line2.substr(43, 8).c_str(), 0)*DE2RA;
    xno = std::strtod(line2.substr(52, 11).c_str(), 0)*(TWOPI / XMNPDA);
    rev_number = std::strtol(line2.substr(63, 5).c_str(), 0, 10);

    // Recover original mean motion (xnodp) and semimajor axis (aodp)
    // from input elements

    double a1 = std::pow(XKE / xno, TOTHRD);
    cosio = std::cos(xincl);
    double theta2 = cosio*cosio;
    x3thm1 = 3.*theta2 - 1.;
    double eosq = eo*eo;
    double betao2 = 1. - eosq;
    double betao = std::sqrt(betao2);
    double del1 = 1.5*CK2*x3thm1 / (a1*a1*betao*betao2);
    double ao = a1*(1. - del1*(.5*TOTHRD + del1*(1. + 134. / 81.*del1)));
    double delo = 1.5*CK2*x3thm1 / (ao*ao*betao*betao2);
    xnodp = xno / (1. + delo);
    aodp = ao / (1. - delo);

    // Initialization

    // For perigee less than 220 kilometers, the isimp flag is set and
    // the equations are truncated to linear variation in sqrt a and
    // quadratic variation in mean anomaly.  Also, the c3 term, the
    // delta omega term, and the delta m term are dropped.

    isimp = ((aodp*(1. - eo) / AE) < (220. / XKMPER + AE));

    // For perigee below 156 km, the values of
    // s and qoms2t are altered

    double s4 = S;
    double qoms24 = QOMS2T;
    double perige = (aodp*(1. - eo) - AE)*XKMPER;
    if (perige < 156.)
    {
        s4 = perige - 78.;
        if (perige <= 98.)
        {
            s4 = 20.;
        }
        qoms24 = std::pow((120. - s4)*AE / XKMPER, 4.);
        s4 = s4 / XKMPER + AE;
    }
    double pinvsq = 1. / (aodp*aodp*betao2*betao2);
    double tsi = 1. / (aodp - s4);
    eta = aodp*eo*tsi;
    double etasq = eta*eta;
    double eeta = eo*eta;
    double psisq = std::abs(1. - etasq);
    double coef = qoms24*std::pow(tsi, 4.);
    double coef1 = coef / std::pow(psisq, 3.5);
    double c2 = coef1*xnodp*(aodp*(1. + 1.5*etasq + eeta*(4. + etasq)) + .75 *
        CK2*tsi / psisq*x3thm1*(8. + 3.*etasq*(8. + etasq)));
    c1 = bstar*c2;
    sinio = std::sin(xincl);
    double a3ovk2 = -XJ3 / CK2*std::pow(AE, 3.);
    double c3 = coef*tsi*a3ovk2*xnodp*AE*sinio / eo;
    x1mth2 = 1. - theta2;
    c4 = 2.*xnodp*coef1*aodp*betao2*(eta *
        (2. + .5*etasq) + eo*(.5 + 2.*etasq) - 2.*CK2*tsi /
        (aodp*psisq)*(-3.*x3thm1*(1. - 2.*eeta + etasq *
        (1.5 - .5*eeta)) + .75*x1mth2*(2.*etasq - eeta *
        (1. + etasq))*std::cos(2.*omegao)));
    c5 = 2.*coef1*aodp*betao2*(1. + 2.75*(etasq + eeta) + eeta*etasq);
    double theta4 = theta2*theta2;
    double temp1 = 3.*CK2*pinvsq*xnodp;
    double temp2 = temp1*CK2*pinvsq;
    double temp3 = 1.25*CK4*pinvsq*pinvsq*xnodp;
    xmdot = xnodp + .5*temp1*betao*x3thm1 + .0625*temp2*betao *
        (13. - 78.*theta2 + 137.*theta4);
    double x1m5th = 1. - 5.*theta2;
    omgdot = -.5*temp1*x1m5th + .0625*temp2*(7. - 114.*theta2 +
        395.*theta4) + temp3*(3. - 36.*theta2 + 49.*theta4);
    double xhdot1 = -temp1*cosio;
    xnodot = xhdot1 + (.5*temp2*(4. - 19.*theta2) + 2.*temp3*(3. -
        7.*theta2))*cosio;
    omgcof = bstar*c3*std::cos(omegao);
    xmcof = -TOTHRD*coef*bstar*AE / eeta;
    xnodcf = 3.5*betao2*xhdot1*c1;
    t2cof = 1.5*c1;
    xlcof = .125*a3ovk2*sinio*(3. + 5.*cosio) / (1. + cosio);
    aycof = .25*a3ovk2*sinio;
    delmo = std::pow(1. + eta*std::cos(xmo), 3.);
    sinmo = std::sin(xmo);
    x7thm1 = 7.*theta2 - 1.;
    if (!isimp)
    {
        double c1sq = c1*c1;
        d2 = 4.*aodp*tsi*c1sq;
        double temp = d2*tsi*c1 / 3.;
        d3 = (17.*aodp + s4)*temp;
        d4 = .5*temp*aodp*tsi*(221.*aodp + 31.*s4)*c1;
        t3cof = d2 + 2.*c1sq;
        t4cof = .25*(3.*d3 + c1*(12.*d2 + 10.*c1sq));
        t5cof = .2*(3.*d4 + 12.*c1*d3 + 6.*d2*d2 + 15.*c1sq*(
            2.*d2 + c1sq));
    }
}

long sgp4::Satellite::get_catalog_number() const
{
    return catalog_number;
}

double sgp4::Satellite::get_epoch() const
{
    return epoch;
}

void sgp4::Satellite::coast(double tsince, double *eci) const
{
    tsince /= 60;

    // Update for secular gravity and atmospheric drag

    double xmdf = xmo + xmdot*tsince;
    double omgadf = omegao + omgdot*tsince;
    double xnoddf = xnodeo + xnodot*tsince;
    double omega = omgadf;
    double xmp = xmdf;
    double tsq = tsince*tsince;
    double xnode = xnoddf + xnodcf*tsq;
    double tempa = 1. - c1*tsince;
    double tempe = bstar*c4*tsince;
    double templ = t2cof*tsq;
    if (!isimp)
    {
        double delomg = omgcof*tsince;
        double delm = xmcof*(std::pow(1. + eta*std::cos(xmdf), 3.) - delmo);
        double temp = delomg + delm;
        xmp = xmdf + temp;
        omega = omgadf - temp;
        double tcube = tsq*tsince;
        double tfour = tsince*tcube;
        tempa = tempa - d2*tsq - d3*tcube - d4*tfour;
        tempe = tempe + bstar*c5*(std::sin(xmp) - sinmo);
        templ = templ + t3cof*tcube +
            tfour*(t4cof + tsince*t5cof);
    }
    double a = aodp*std::pow(tempa, 2.);
    double e = eo - tempe;
    double xl = xmp + omega + xnode + xnodp*templ;
    double beta = std::sqrt(1. - e*e);
    double xn = XKE / std::pow(a, 1.5);

    // Long period periodics

    double axn = e*std::cos(omega);
    double temp = 1. / (a*beta*beta);
    double xll = temp*xlcof*axn;
    double aynl = temp*aycof;
    double xlt = xl + xll;
    double ayn = e*std::sin(omega) + aynl;

    // Solve Keplers equation

    double capu = std::fmod(xlt - xnode, TWOPI);
    if (capu < 0)
    {
        capu += TWOPI;
    }
    double temp2 = capu,
        temp3 = 0, temp4 = 0, temp5 = 0, temp6 = 0,
        sinepw = 0, cosepw = 0;
    for (int i = 0; i < 10; ++i)
    {
        sinepw = std::sin(temp2);
        cosepw = std::cos(temp2);
        temp3 = axn*sinepw;
        temp4 = ayn*cosepw;
        temp5 = axn*cosepw;
        temp6 = ayn*sinepw;
        double epw =
            (capu - temp4 + temp3 - temp2) / (1. - temp5 - temp6) + temp2;
        if (std::abs(epw - temp2) <= E6A)
        {
            break;
        }
        temp2 = epw;
    }

    // Short period preliminary quantities

    double ecose = temp5 + temp6;
    double esine = temp3 - temp4;
    double elsq = axn*axn + ayn*ayn;
    temp = 1. - elsq;
    double pl = a*temp;
    double r = a*(1. - ecose);
    double temp1 = 1. / r;
    double rdot = XKE*std::sqrt(a)*esine*temp1;
    double rfdot = XKE*std::sqrt(pl)*temp1;
    temp2 = a*temp1;
    double betal = std::sqrt(temp);
    temp3 = 1. / (1. + betal);
    double cosu = temp2*(cosepw - axn + ayn*esine*temp3);
    double sinu = temp2*(sinepw - ayn - axn*esine*temp3);
    double u = std::atan2(sinu, cosu);
    double sin2u = 2.*sinu*cosu;
    double cos2u = 2.*cosu*cosu - 1.;
    temp = 1. / pl;
    temp1 = CK2*temp;
    temp2 = temp1*temp;

    // Update for short periodics

    double rk = r*(1. - 1.5*temp2*betal*x3thm1) + .5*temp1*x1mth2*cos2u;
    double uk = u - .25*temp2*x7thm1*sin2u;
    double xnodek = xnode + 1.5*temp2*cosio*sin2u;
    double xinck = xincl + 1.5*temp2*cosio*sinio*cos2u;
    double rdotk = rdot - xn*temp1*x1mth2*sin2u;
    double rfdotk = rfdot + xn*temp1*(x1mth2*cos2u + 1.5*x3thm1);

    // Orientation vectors

    double sinuk = std::sin(uk);
    double cosuk = std::cos(uk);
    double sinik = std::sin(xinck);
    double cosik = std::cos(xinck);
    double sinnok = std::sin(xnodek);
    double cosnok = std::cos(xnodek);
    double xmx = -sinnok*cosik;
    double xmy = cosnok*cosik;
    double ux = xmx*sinuk + cosnok*cosuk;
    double uy = xmy*sinuk + sinnok*cosuk;
    double uz = sinik*sinuk;
    double vx = xmx*cosuk - cosnok*sinuk;
    double vy = xmy*cosuk - sinnok*sinuk;
    double vz = sinik*cosuk;

    // Position and velocity

    eci[0] = (rk*ux)*POS_UNITS;
    eci[1] = (rk*uy)*POS_UNITS;
    eci[2] = (rk*uz)*POS_UNITS;
    eci[3] = (rdotk*ux + rfdotk*vx)*VEL_UNITS;
    eci[4] = (rdotk*uy + rfdotk*vy)*VEL_UNITS;
    eci[5] = (rdotk*uz + rfdotk*vz)*VEL_UNITS;
}

void sgp4::Satellite::eci_to_efg(
    double jd, const double *eci, double *efg) const
{
    // Greenwich Mean Sidereal Time (radians)
    double t = (jd - 2451545) / 36525;
    double gmst = std::fmod((67310.54841 +
        t * (3164400184.812866 + t * (0.093104 - t * 6.2e-6))) * TWOPI / 86400,
        TWOPI);

    double s = std::sin(gmst);
    double c = std::cos(gmst);
    efg[0] =  c * eci[0] + s * eci[1];
    efg[1] = -s * eci[0] + c * eci[1];
    efg[2] = eci[2];
    efg[3] =  OMEGA * efg[1] + c * eci[3] + s * eci[4];
    efg[4] = -OMEGA * efg[0] - s * eci[3] + c * eci[4];
    efg[5] = eci[5];
}
