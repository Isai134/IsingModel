using Random
using Statistics
using Plots


function initialstate(N)
    state = 2 .* rand(Bool, (N, N)) .- 1
    return state
end

function mcmove(config, beta)
    N = size(config, 1)
    for _ in 1:N, _ in 1:N
        a, b = rand(1:N), rand(1:N)
        a, b = Int(a), Int(b)
        s = config[a, b]
        nb = config[mod1(a+1, N), b] + config[a, mod1(b+1, N)] + config[mod1(a-1, N), b] + config[a, mod1(b-1, N)]
        cost = 2s * nb
        if cost < 0 || rand() < exp(-cost * beta)
            s *= -1
        end
        config[a, b] = s
    end
    return config
end

function calcEnergy(config)
    energy = 0
    N = size(config, 1)
    for i in 1:N, j in 1:N
        S = config[i, j]
        nb = config[mod1(i+1, N), j] + config[i, mod1(j+1, N)] + config[mod1(i-1, N), j] + config[i, mod1(j-1, N)]
        energy += -nb * S
    end
    return energy / 4
end

function calcMag(config)
    return sum(config)
end
nt = 88 
N = 8
eqSteps = 200000
mcSteps = 200000

T = range(1.53, stop=3.28, length=nt)
E = zeros(nt)
M = zeros(nt)
C = zeros(nt)
X = zeros(nt)
n1, n2 = 1.0 / (mcSteps * N * N), 1.0 / (mcSteps * mcSteps * N * N)

for tt in 1:nt
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N)
    iT = 1.0 / T[tt]
    iT2 = iT * iT

    for i in 1:eqSteps
        mcmove(config, iT)
    end

    for i in 1:mcSteps
        mcmove(config, iT)
        Ene = calcEnergy(config)
        Mag = calcMag(config)

        E1 += Ene
        M1 += Mag
        M2 += Mag * Mag
        E2 += Ene * Ene
    end

    E[tt] = n1 * E1
    M[tt] = n1 * M1
    C[tt] = (n1 * E2 - n2 * E1 * E1) * iT2
    X[tt] = (n1 * M2 - n2 * M1 * M1) * iT
end
Nf = 8
eqStepsf = 200000
mcStepsf = 200000

Tf = range(1.53, stop=3.28, length=nt)
Ef = zeros(nt)
Mf = zeros(nt)
Cf = zeros(nt)
Xf = zeros(nt)
n1, n2 = 1.0 / (mcStepsf * Nf * Nf), 1.0 / (mcStepsf * mcStepsf * Nf * Nf)

for tt in 1:nt
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N)
    iT = 1.0 / T[tt]
    iT2 = iT * iT

    for i in 1:eqStepsf
        mcmove(config, iT)
    end

    for i in 1:mcStepsf
        mcmove(config, iT)
        Ene = calcEnergy(config)
        Mag = calcMag(config)

        E1 += Ene
        M1 += Mag
        M2 += Mag * Mag
        E2 += Ene * Ene
    end

    Ef[tt] = n1 * E1
    Mf[tt] = n1 * M1
    Cf[tt] = (n1 * E2 - n2 * E1 * E1) * iT2
    Xf[tt] = (n1 * M2 - n2 * M1 * M1) * iT
end

scatter(T, E, markersize= 4, xlabel="Temperatura (T)", ylabel="Energía",tickfontsize=8,guidefontsize=12, alpha = 0.4, label="Pasos=200");
scatter!(Tf, Ef, markersize= 4, color=:IndianRed, xlabel="Temperatura (T)", ylabel="Energía",tickfontsize=8,guidefontsize=12, label="Pasos=200k")

scatter(T, abs.(M), markersize=4, xlabel="Temperatura (T)", ylabel="Magnetización",tickfontsize=8,guidefontsize=12, alpha = 0.4, label="Pasos=200");
scatter!(Tf, abs.(Mf), markersize=4, color=:RoyalBlue, xlabel="Temperatura (T)", ylabel="Magnetización",tickfontsize=8,guidefontsize=12, label="Pasos=200k")

scatter(T, C, markersize=4, color=:IndianRed, xlabel="Temperatura (T)", ylabel="Calor específico",tickfontsize=8,guidefontsize=12, alpha = 0.4, label="Pasos=200");
scatter!(Tf, Cf, markersize=4, color=:IndianRed, xlabel="Temperatura (T)", ylabel="Calor específico",tickfontsize=8,guidefontsize=12, label="Pasos=200k")

scatter(T, X, markersize=4, xlabel="Temperatura (T)", ylabel="Susceptibilidad",tickfontsize=8,guidefontsize=12, alpha = 0.4, label="Pasos=200");
scatter!(Tf, Xf, markersize=4, color=:RoyalBlue, xlabel="Temperatura (T)", ylabel="Susceptibilidad",tickfontsize=8,guidefontsize=12, label="Pasos=200k")

plot(plt_E, plt_M, plt_C, plt_X, layout=(2, 2), size=(1200, 800))


mutable struct Ising
end

function mcmove(config, N, beta)
    for _ in 1:N, _ in 1:N
        a, b = rand(1:N), rand(1:N)
        a, b = Int(a), Int(b)
        s = config[a, b]
        nb = config[mod1(a+1, N), b] + config[a, mod1(b+1, N)] + config[mod1(a-1, N), b] + config[a, mod1(b-1, N)]
        cost = 2s * nb
        if cost < 0 || rand() < exp(-cost * beta)
            s *= -1
        end
        config[a, b] = s
    end
    return config
end

function simulate()
    N, temp = 64, 0.4
    config = 2 .* rand(Bool, (N, N)) .- 1
    f = plot(layout=(3, 2), size=(900, 900))
    configPlot(f, config, 0, N, 1)

    msrmnt = 1001
    for i in 0:msrmnt-1
        mcmove(config, N, 1.0 / temp)
        if i == 1      configPlot(f, config, i, N, 2)
        elseif i == 4  configPlot(f, config, i, N, 3)
        elseif i == 32 configPlot(f, config, i, N, 4)
        elseif i == 100 configPlot(f, config, i, N, 5)
        elseif i == 1000 configPlot(f, config, i, N, 6)
        end
    end
    return f
end

function configPlot(f, config, i, N, n_)
    x = 1:N
    y = 1:N
    heatmap!(x, y, config, subplot=n_, c=:RdBu, title="Tiempo = $i s", color=:auto, yticks=false, xticks=false)
end

rm = Ising()
simulate()
