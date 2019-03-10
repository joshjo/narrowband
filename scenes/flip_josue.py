from manta import *

dim = 2
particle_number = 2
res = 32
gs = vec3(res, res, res)

narrowBandWidth = 3
combineBandWidth = narrowBandWidth - 1

if dim == 2:
    gs.z = 1
    particleNumber = 3

s = Solver(name='Josue', gridSize=gs, dim=dim)
s.timestep = 0.5

minParticles = pow(2, dim)

flags   = s.create(FlagGrid)

phi     = s.create(LevelsetGrid)
phiParts = s.create(LevelsetGrid)

vel     = s.create(MACGrid)
velOld  = s.create(MACGrid)
velParts = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3 = s.create(MACGrid)
pp      = s.create(BasicParticleSystem)

pVel    = pp.create(PdataVec3)
mesh      = s.create(Mesh)

pindex = s.create(ParticleIndexSystem)
gpi    = s.create(IntGrid)

flags.initDomain(boundaryWidth=0)
phi.initFromFlags(flags)

fluidbox = Box(parent=s, p0=gs*vec3(0, 0, 0), p1=gs*vec3(0.2, 1, 1))
phi.join(fluidbox.computeLevelset())

# apic part
mass    = s.create(MACGrid)
pCx     = pp.create(PdataVec3)
pCy     = pp.create(PdataVec3)
pCz     = pp.create(PdataVec3)

# fluidDam = Box(parent=s, p0=gs*vec3(0, 0.15, 0), p1=gs*vec3(0.4, 0.5, 0.8))
# phi.join(fluidDam.computeLevelset())

flags.updateFromLevelset(phi)

sampleLevelsetWithParticles(
    phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.1)

mapGridToPartsVec3(source=vel, parts=pp, target=pVel)


if GUI:
    gui = Gui()
    gui.show()
    gui.pause()

step = -1

for i in range(2000):
    maxVel = vel.getMax()
    s.adaptTimestep(maxVel)

    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK2, deleteInObstacle=False)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)
    flags.updateFromLevelset(phi)

    # # Here insert the correct particle positions

    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

    gridParticleIndex(parts=pp , flags=flags, indexSys=pindex, index=gpi)
    unionParticleLevelset(pp, pindex, flags, gpi, phiParts)

    phi.addConst(1.)
    phi.join(phiParts)
    extrapolateLsSimple(phi=phi, distance=narrowBandWidth+2, inside=True)

    extrapolateLsSimple(phi=phi, distance=3)
    flags.updateFromLevelset(phi)

    # APIC Transfer to MAC
    particleMACTransfers(flags=flags, vel=velParts, parts=pp, partVel=pVel, cpx=pCx, cpy=pCy, cpz=pCz, mass=mass)
    combineGridVel(vel=velParts, weight=tmpVec3, combineVel=vel, phi=phi, narrowBand=combineBandWidth, thresh=0)
    velOld.copyFrom(vel)

    addGravity(flags=flags, vel=vel, gravity=(0, -0.003, 0))
    setWallBcs(flags=flags, vel=vel)
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setWallBcs(flags=flags, vel=vel)

    extrapolateMACSimple(flags=flags, vel=vel, distance=(int(maxVel*1.25 + 2.)))

    # APIC Transfer to vel
    particleGridTransfers(cpx=pCx, cpy=pCy, cpz=pCz, partVel=pVel, parts=pp, vel=vel, flags=flags)

    if dim==3:
        phi.createMesh(mesh)

    pVel.setSource(vel, isMAC=True)
    phi.setBoundNeumann(0)
    adjustNumber(
        parts=pp, vel=vel, flags=flags,
        minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi,
        narrowBand=narrowBandWidth
    )

    s.step()
