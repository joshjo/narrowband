#include "particle.h"
#include "grid.h"


namespace Manta {


KERNEL(pts, single)
void particleMACTransfersKernel(
	const BasicParticleSystem &p
	, MACGrid &mg
	, MACGrid &vg
	, const ParticleDataImpl<Vec3> &vp
	, const ParticleDataImpl<Vec3> &cpx
	, const ParticleDataImpl<Vec3> &cpy
	, const ParticleDataImpl<Vec3> &cpz
	, const ParticleDataImpl<int> *ptype
	, const int exclude
) {
	if (!p.isActive(idx) || (ptype &&	 ((*ptype)[idx] & exclude))) {
		return;
	}
	IndexInt dX[2] = { 0, vg.getStrideX() };
	IndexInt dY[2] = { 0, vg.getStrideY() };
	IndexInt dZ[2] = { 0, vg.getStrideZ() };

	const Vec3 &pos = p[idx].pos;
	const Vec3 &vel = vp[idx];

	IndexInt fi = (IndexInt)pos.x;
	IndexInt fj = (IndexInt)pos.y;
	IndexInt fk = (IndexInt)pos.z;

	IndexInt ci = (IndexInt)(pos.x - 0.5);
	IndexInt cj = (IndexInt)(pos.y - 0.5);
	IndexInt ck = (IndexInt)(pos.z - 0.5);

    Real wfi = pos.x - fi;
    Real wfj = pos.y - fj;
    Real wfk = pos.z - fk;

    Real wci = pos.x - ci - 0.5;
    Real wcj = pos.y - cj - 0.5;
    Real wck = pos.z - ck - 0.5;

	IndexInt fgidx = fi*dX[1] + cj*dY[1] + ck*dZ[1];
	Vec3 fgpos(fi, cj+0.5, ck+0.5);
	Real iwi[2] = { 1 - wfi, wfi };
	Real iwj[2] = { 1 - wcj, wcj };
	Real iwk[2] = { 1 - wck, wck };
	for (int i=0; i<2; ++i) {
		for (int j=0; j<2; ++j) {
			for (int k=0; k<2; ++k) {
				const Real w = iwi[i] * iwj[j] * iwk[k];
				mg[fgidx + dX[i] + dY[j] + dZ[k]].x += w;
				vg[fgidx + dX[i] + dY[j] + dZ[k]].x += w * vel.x;
				vg[fgidx + dX[i] + dY[j] + dZ[k]].x += w * dot(cpx[idx], fgpos + Vec3(i, j, k) - pos);
			}
		}
	}
	IndexInt cgidx = ci * dX[1] + fj*dY[1] + ck*dZ[1];
	Vec3 cgpos(ci + 0.5, fj, ck+0.5);
	Real jwi[2] = { 1 - wci, wci };
	Real jwj[2] = { 1 - wfj, wfj };
	Real jwk[2] = { 1 - wck, wck };
	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < 2; ++j) {
			for(int k = 0; k < 2; ++k) {
				const Real w = jwi[i] * jwj[j] * jwk[k];
				mg[cgidx + dX[i] + dY[j] + dZ[k]].y += w;
				vg[cgidx + dX[i] + dY[j] + dZ[k]].y += w * vel.y;
				vg[cgidx + dX[i] + dY[j] + dZ[k]].y += w * dot(cpy[idx], cgpos + Vec3(i, j, k) - pos);
			}
		}
	}

	if(vg.is3D()) {
		const IndexInt gidx = ci*dX[1] + cj*dY[1] + fk*dZ[1];
		const Vec3 gpos(ci+0.5, cj+0.5, fk);
		const Real wi[2] = { Real(1)-wci, wci };
		const Real wj[2] = { Real(1)-wcj, wcj };
		const Real wk[2] = { Real(1)-wfk, wfk };
		for(int i=0; i<2; ++i)
			for(int j=0; j<2; ++j)
				for(int k=0; k<2; ++k) {
					const Real w = wi[i]*wj[j]*wk[k];
					mg[gidx+dX[i]+dY[j]+dZ[k]].z += w;
					vg[gidx+dX[i]+dY[j]+dZ[k]].z += w*vel.z;
					vg[gidx+dX[i]+dY[j]+dZ[k]].z += w*dot(cpz[idx], gpos + Vec3(i, j, k) - pos);
				}
	}
}

PYTHON()
void particleMACTransfers(
	const FlagGrid &flags
	, MACGrid &vel
	, const BasicParticleSystem &parts
	, const ParticleDataImpl<Vec3> &partVel
	, const ParticleDataImpl<Vec3> &cpx
	, const ParticleDataImpl<Vec3> &cpy
	, const ParticleDataImpl<Vec3> &cpz
	, MACGrid *mass=NULL
	, const ParticleDataImpl<int> *ptype=NULL
	, const int exclude=0
) {
	const bool freeMass = !mass;
	if(!mass) {
		mass = new MACGrid(flags.getParent());
	} else {
		mass->clear();
	}

	vel.clear();
	particleMACTransfersKernel(parts, *mass, vel, partVel, cpx, cpy, cpz, ptype, exclude);
	mass->stomp(VECTOR_EPSILON);
	vel.safeDivide(*mass);

	if (freeMass) {
		delete mass;
	}
}

KERNEL(pts)
void particleGridTransfersKernel(
	ParticleDataImpl<Vec3> &vp
    , ParticleDataImpl<Vec3> &cpx
    , ParticleDataImpl<Vec3> &cpy
    , ParticleDataImpl<Vec3> &cpz
    , const BasicParticleSystem &p
    , const MACGrid &vg, const FlagGrid &flags
    , const ParticleDataImpl<int> *ptype
    , const int exclude
)
{
	if (!p.isActive(idx) || (ptype && ((*ptype)[idx] & exclude))) {
		return;
	}
	vp[idx] = cpx[idx] = cpy[idx] = cpz[idx] = Vec3(Real(0));
	const IndexInt dX[2] = {0, vg.getStrideX()}, dY[2] = {0, vg.getStrideY()}, dZ[2] = {0, vg.getStrideZ()};
	const Real gw[2] = {-Real(1), Real(1)};

	const Vec3 &pos = p[idx].pos;

	IndexInt fi = (IndexInt)pos.x;
	IndexInt fj = (IndexInt)pos.y;
	IndexInt fk = (IndexInt)pos.z;

	IndexInt ci = (IndexInt)(pos.x - 0.5);
	IndexInt cj = (IndexInt)(pos.y - 0.5);
	IndexInt ck = (IndexInt)(pos.z - 0.5);

    Real wfi = pos.x - fi;
    Real wfj = pos.y - fj;
    Real wfk = pos.z - fk;

    Real wci = pos.x - ci - 0.5;
    Real wcj = pos.y - cj - 0.5;
    Real wck = pos.z - ck - 0.5;

	IndexInt gidx = fi * dX[1] + cj * dY[1] + ck * dZ[1];
	Real iwx[2] = { Real(1) - wfi, wfi };
	Real iwy[2] = { Real(1) - wcj, wcj };
	Real iwz[2] = { Real(1) - wck, wck };
	for(int i = 0; i < 2; ++i) {
		for(int j = 0; j < 2; ++j) {
			for(int k = 0; k < 2; ++k) {
				const IndexInt vidx = gidx + dX[i] + dY[j] + dZ[k];
				Real vgx = vg[vidx].x;
				vp[idx].x  += iwx[i] * iwy[j] * iwz[k] * vgx;
				cpx[idx].x +=  gw[i] * iwy[j] * iwz[k] * vgx;
				cpx[idx].y += iwx[i] *  gw[j] * iwz[k] * vgx;
				cpx[idx].z += iwx[i] * iwy[j] *  gw[k] * vgx;
			}
		}
	}

	gidx = ci * dX[1] + fj * dY[1] + ck * dZ[1];
	Real jwx[2] = { Real(1)-wci, wci };
	Real jwy[2] = { Real(1)-wfj, wfj };
	Real jwz[2] = { Real(1)-wck, wck };
	for(int i=0; i < 2; ++i) {
		for(int j=0; j < 2; ++j) {
			for(int k=0; k < 2; ++k) {
				const IndexInt vidx = gidx + dX[i] + dY[j] + dZ[k];
				Real vgy = vg[vidx].y;
				vp[idx].y  += jwx[i] * jwy[j] * jwz[k] * vgy;
				cpy[idx].x +=  gw[i] * jwy[j] * jwz[k] * vgy;
				cpy[idx].y += jwx[i] *  gw[j] * jwz[k] * vgy;
				cpy[idx].z += jwx[i] * jwy[j] *  gw[k] * vgy;
			}
		}
	}
	if(vg.is3D()) {
		const IndexInt gidx = ci * dX[1] + cj * dY[1] + fk * dZ[1];
		const Real kwx[2] = { Real(1) - wci, wci };
		const Real kwy[2] = { Real(1) - wcj, wcj };
		const Real kwz[2] = { Real(1) - wfk, wfk };
		for(int i=0; i<2; ++i)
			for(int j=0; j<2; ++j)
				for(int k=0; k<2; ++k) {
					const IndexInt vidx = gidx + dX[i] + dY[j] + dZ[k];
					Real vgz = vg[vidx].z;
					vp[idx].z  += kwx[i] * kwy[j] * kwz[k] * vgz;
					cpz[idx].x += gw[i] * kwy[j] * kwz[k] * vgz;
					cpz[idx].y += kwx[i] * gw[j] * kwz[k] * vgz;
					cpz[idx].z += kwx[i] * kwy[j] * gw[k] * vgz;
				}
	}
}

PYTHON()
void particleGridTransfers(
    ParticleDataImpl<Vec3> &cpx
    , ParticleDataImpl<Vec3> &cpy
    , ParticleDataImpl<Vec3> &cpz
    , ParticleDataImpl<Vec3> &partVel
    , const BasicParticleSystem &parts
    , const MACGrid &vel
    , const FlagGrid &flags
    , const ParticleDataImpl<int> *ptype=NULL
	, const int exclude=0
) {

    particleGridTransfersKernel(cpx, cpy, cpz, partVel, parts, vel, flags, ptype, exclude);
}

}
