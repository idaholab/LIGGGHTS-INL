/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Andreas Aigner (DCS Computing GmbH, Linz)

    Copyright 2018-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_VERSION_STRUCT_H
#define LMP_VERSION_STRUCT_H

struct Version
{
    int major_;
    int minor_;

    Version(int major=0, int minor=0)
        : major_(major), minor_(minor)
    {

    }

    Version& operator=(const Version& a)
    {
        major_ = a.major_;
        minor_ = a.minor_;
        return *this;
    }

    bool operator==(const Version& a) const
    {
        return (major_ == a.major_ && minor_ == a.minor_);
    }

    bool operator!=(const Version& a) const
    {
        return (major_ != a.major_ || minor_ != a.minor_);
    }

    bool operator<(const Version& a) const
    {
        return (major_ < a.major_ || (major_ == a.major_ && minor_ < a.minor_));
    }

    bool operator>(const Version& a) const
    {
        return (major_ > a.major_ || (major_ == a.major_ && minor_ > a.minor_));
    }

    bool operator<=(const Version& a) const
    {
        return *this<a || *this==a;
    }

    bool operator>=(const Version& a) const
    {
        return *this>a || *this==a;
    }
};

#endif
