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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_ASSOCIATIVE_POINTER_ARRAY_H
#define LMP_ASSOCIATIVE_POINTER_ARRAY_H

#include <string>
#include <list>
#include <algorithm>
#include "memory.h"

namespace LAMMPS_NS
{
  #define ID_LEN 100

template<typename T>
class AssociativePointerArray
{
      public:
        AssociativePointerArray();
        ~AssociativePointerArray();

        template <typename U>
        U* add(const char *_id, const char* _comm, const char* _ref, const char *_restart,int _scalePower = 1);

        void remove(const char *_id);

        template <typename U>
        U* getPointerById(const char *_id);

        T* getBasePointerById(const char *_id);

        template <typename U>
        U* getPointerByIndex(int i);

        T* getBasePointerByIndex(int i) const;

        int size() const;

        bool sameLength(int _len);

        inline void copyElement(int from, int to);
        inline void addUninitializedElement();
        inline void addZeroElement();
        inline void deleteAllElements();
        inline void deleteRestart(bool scale,bool translate,bool rotate);
        inline void deleteElement(int n);
        inline void deleteForwardElement(int n,bool scale,bool translate,bool rotate);
        inline void deleteRestartElement(int n,bool scale,bool translate,bool rotate);
        inline void deleteRestartGlobal(bool scale,bool translate,bool rotate);

        inline void clearReverse(bool scale,bool translate,bool rotate);

        inline bool calcStatistics();
        inline int  maxStatLevel() const;

        inline void storeOrig(class AssociativePointerArray &orig);
        inline void storeOrig(const char *_id,class AssociativePointerArray &orig);
        inline bool reset(class AssociativePointerArray &orig);
        inline bool reset(const char *_id,class AssociativePointerArray &orig);

        void rotate(const double * const dQ);
        void move(const double * const delta);
        void moveElement(const int i, const double * const delta);
        void scale(double factor);

        inline int bufSize(int operation,bool scale,bool translate,bool rotate) const;
        inline int pushToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate) const;
        inline int popFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate) const;
        inline int pushElemListToBuffer(const int n, const int *const list, double *const buf, const int operation, std::list<std::string> *const properties, const int pbc_flag, const int *const pbc, Domain *const domain, const bool scale, const bool translate, const bool rotate) const;
        inline int popElemListFromBuffer(int first, int n, double *buf, int operation, std::list<std::string> * properties, bool scale,bool translate, bool rotate);
        inline int pushElemListToBufferReverse(int first, int n, double *buf, int operation, std::list<std::string> *properties,bool scale,bool translate, bool rotate) const;
        inline int popElemListFromBufferReverse(int n, const int *const list, double *buf, int operation, std::list<std::string> *properties, bool scale,bool translate, bool rotate);

        inline int elemBufSize(int operation, std::list<std::string> * properties, bool scale,bool translate,bool rotate) const;
        inline int pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate) const;
        inline int popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        int idToIndex(const char *_id);
        void indexToId(int index, char *_id);

      private:

        T **content_;
        int numElem_, maxElem_;

        void growArrays();
};

  // *************************************
  #include "associative_pointer_array_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* ASSOCIATIVEPOINTERARRAY_H_ */
