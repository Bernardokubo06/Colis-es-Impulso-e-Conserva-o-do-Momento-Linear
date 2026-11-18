# ADO-3 - Colisões, impulso e momento linear (Opção A: colisões elásticas)
# Aluno(s): ______________________________________

import itertools

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, RadioButtons, Slider

X1_0, X2_0 = -6.0, 2.0
XMIN, XMAX = -10.0, 10.0
T_MAX = 20.0
LADO_POR_RAIZ = 0.8

INTERVALO_MS = 30
ESCALA_TEMPO = 0.5   # câmera lenta
PASSO_T = INTERVALO_MS / 1000.0 * ESCALA_TEMPO

# (m1, m2, v1), com o alvo em repouso
CASOS = {
    "Caso 1: m1 = m2": (1.0, 1.0, 4.0),
    "Caso 2: m1 < m2": (1.0, 4.0, 5.0),
    "Caso 3: m1 > m2": (4.0, 1.0, 5.0),
}

COR1, COR2 = "#2a6fdb", "#e8832a"
COR_ANTES, COR_DEPOIS = "#8d99ae", "#2a9d8f"


# ---------------------------------------------------------------- física
def momento(m1, v1, m2, v2):
    return m1 * v1 + m2 * v2


def cinetica(m1, v1, m2, v2):
    return 0.5 * m1 * v1 ** 2 + 0.5 * m2 * v2 ** 2


def colisao_elastica(m1, m2, v1, v2):
    soma = m1 + m2
    u1 = ((m1 - m2) * v1 + 2.0 * m2 * v2) / soma
    u2 = ((m2 - m1) * v2 + 2.0 * m1 * v1) / soma
    return u1, u2


def impulso(m, v, u):
    return m * (u - v)


def forca_media(J, dt):
    return abs(J) / dt


def lado_bloco(m):
    return LADO_POR_RAIZ * np.sqrt(m)


def instante_choque(x1, x2, meia_soma, v1, v2):
    # só há choque se o bloco 1 se aproxima do 2
    if v1 <= v2:
        return None
    return (x2 - x1 - meia_soma) / (v1 - v2)


def _tempo_de_saida(x0, v, lado, t0=0.0):
    if v > 0:
        return t0 + (XMAX - lado / 2 - x0) / v
    if v < 0:
        return t0 + (XMIN + lado / 2 - x0) / v
    return np.inf


def calcular_cenario(m1, m2, v1, v2, dt_ms):
    dt = dt_ms * 1e-3
    w1, w2 = lado_bloco(m1), lado_bloco(m2)
    tc = instante_choque(X1_0, X2_0, (w1 + w2) / 2, v1, v2)
    colide = tc is not None

    if colide:
        u1, u2 = colisao_elastica(m1, m2, v1, v2)
    else:
        u1, u2 = v1, v2

    J1, J2 = impulso(m1, v1, u1), impulso(m2, v2, u2)
    p_a, p_d = momento(m1, v1, m2, v2), momento(m1, u1, m2, u2)
    K_a, K_d = cinetica(m1, v1, m2, v2), cinetica(m1, u1, m2, u2)

    x1c = X1_0 + v1 * tc if colide else None
    x2c = X2_0 + v2 * tc if colide else None

    # a simulação termina quando um bloco sai da janela ou em T_MAX
    ex1, ex2 = _tempo_de_saida(X1_0, v1, w1), _tempo_de_saida(X2_0, v2, w2)
    fora = False
    if colide:
        if tc < min(ex1, ex2, T_MAX):
            ex1 = _tempo_de_saida(x1c, u1, w1, tc)
            ex2 = _tempo_de_saida(x2c, u2, w2, tc)
        else:
            fora = True
    t_fim = min(T_MAX, ex1, ex2)

    return dict(m1=m1, m2=m2, v1=v1, v2=v2, u1=u1, u2=u2, dt_ms=dt_ms,
                w1=w1, w2=w2, tc=tc, colide=colide, fora=fora, t_fim=t_fim,
                x1c=x1c, x2c=x2c, J1=J1, J2=J2,
                F_med=forca_media(J1, dt),
                p_antes=p_a, p_depois=p_d,
                K_antes=K_a, K_depois=K_d, dK=K_d - K_a)


def estado_em(c, t):
    if c["colide"] and t >= c["tc"]:
        tau = t - c["tc"]
        return (c["x1c"] + c["u1"] * tau, c["x2c"] + c["u2"] * tau,
                c["u1"], c["u2"])
    return X1_0 + c["v1"] * t, X2_0 + c["v2"] * t, c["v1"], c["v2"]


# ------------------------------------------------------------- interface
def formata(x, casas=2, sinal=True):
    x = 0.0 if abs(x) < 1e-9 else x  # evita "-0.00"
    return f"{x:+.{casas}f}" if sinal else f"{x:.{casas}f}"


class Interface:
    FIG_W, FIG_H = 13.0, 8.6
    AX_ANIM = [0.05, 0.60, 0.57, 0.33]

    def __init__(self):
        self.t = 0.0
        self.pausado = False
        self.cen = None
        self._ajustando = False
        self._montar_figura()
        self._montar_controles()
        self.reiniciar()
        self.anim = FuncAnimation(self.fig, self.atualizar,
                                  frames=itertools.count(),
                                  interval=INTERVALO_MS, blit=False,
                                  cache_frame_data=False)

    def _montar_figura(self):
        self.fig = plt.figure(figsize=(self.FIG_W, self.FIG_H))
        try:
            self.fig.canvas.manager.set_window_title(
                "ADO-3 - Colisões elásticas (Opção A)")
        except Exception:
            pass
        self.fig.suptitle("ADO-3 - Colisão frontal elástica: impulso e "
                          "conservação do momento linear",
                          fontsize=13, fontweight="bold", y=0.985)

        ax = self.ax = self.fig.add_axes(self.AX_ANIM)
        razao = (self.AX_ANIM[3] * self.FIG_H) / (self.AX_ANIM[2] * self.FIG_W)
        ax.set_xlim(XMIN, XMAX)
        ax.set_ylim(-1.0, -1.0 + (XMAX - XMIN) * razao)  # escala 1:1
        ax.set_yticks([])
        ax.set_xlabel("posição x (m)")
        ax.set_facecolor("#f7f9fc")
        ax.axhline(0, color="k", lw=2)

        self.r1 = Rectangle((0, 0), 1, 1, fc=COR1, ec="k", lw=1.2)
        self.r2 = Rectangle((0, 0), 1, 1, fc=COR2, ec="k", lw=1.2)
        ax.add_patch(self.r1)
        ax.add_patch(self.r2)
        kw = dict(ha="center", va="center", color="white",
                  fontweight="bold", fontsize=10)
        self.n1 = ax.text(0, 0, "1", **kw)
        self.n2 = ax.text(0, 0, "2", **kw)

        kw = dict(transform=ax.transAxes, va="top", fontsize=9,
                  family="monospace")
        self.info1 = ax.text(0.01, 0.97, "", ha="left", color=COR1, **kw)
        self.info2 = ax.text(0.99, 0.97, "", ha="right", color=COR2, **kw)
        self.status = ax.text(0.5, 0.97, "", ha="center", color="k",
                              fontweight="bold", **{**kw, "family": "sans-serif"})
        self.aviso = ax.text(0.5, 0.62, "", transform=ax.transAxes,
                             ha="center", va="center", color="#c0392b",
                             fontsize=12, fontweight="bold",
                             bbox=dict(fc="white", ec="#c0392b", boxstyle="round"))

        self.ax_p, self.bar_p, self.rot_p = self._montar_barras(
            [0.68, 0.71, 0.12, 0.21], "Momento total p (kg·m/s)")
        self.ax_k, self.bar_k, self.rot_k = self._montar_barras(
            [0.86, 0.71, 0.12, 0.21], "Energia cinética total K (J)")

        self.ax_txt = self.fig.add_axes([0.655, 0.04, 0.335, 0.59])
        self.ax_txt.set_xticks([])
        self.ax_txt.set_yticks([])
        self.ax_txt.set_facecolor("#fbfbfb")
        self.painel = self.ax_txt.text(0.03, 0.98, "", va="top", ha="left",
                                       family="monospace", fontsize=10,
                                       transform=self.ax_txt.transAxes)

    def _montar_barras(self, rect, titulo):
        ax = self.fig.add_axes(rect)
        barras = ax.bar([0, 1], [0, 0], width=0.6, edgecolor="k",
                        color=[COR_ANTES, COR_DEPOIS])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["antes", "depois"])
        ax.set_title(titulo, fontsize=9)
        ax.axhline(0, color="k", lw=0.8)
        rotulos = [ax.text(i, 0, "", ha="center", fontsize=9,
                           fontweight="bold") for i in (0, 1)]
        return ax, barras, rotulos

    def _montar_controles(self):
        defs = [
            ("m1 (kg)", 0.5, 10.0, 1.0, 0.1, COR1),
            ("m2 (kg)", 0.5, 10.0, 1.0, 0.1, COR2),
            ("v1 (m/s)", -10.0, 10.0, 4.0, 0.1, COR1),
            ("v2 (m/s)", -10.0, 10.0, 0.0, 0.1, COR2),
            ("Δt contato (ms)", 1.0, 100.0, 10.0, 1.0, "#6c757d"),
        ]
        self.sliders = []
        for i, (nome, vmin, vmax, v0, passo, cor) in enumerate(defs):
            ax = self.fig.add_axes([0.15, 0.51 - 0.055 * i, 0.40, 0.028])
            fmt = "%.0f" if passo >= 1 else "%.1f"
            s = Slider(ax, nome, vmin, vmax, valinit=v0, valstep=passo,
                       valfmt=fmt, color=cor)
            s.on_changed(self._ao_mudar)
            self.sliders.append(s)
        self.s_m1, self.s_m2, self.s_v1, self.s_v2, self.s_dt = self.sliders

        ax_r = self.fig.add_axes([0.05, 0.06, 0.26, 0.15])
        ax_r.set_title("Caso (opção A)", fontsize=10, loc="left")
        self.radio = RadioButtons(ax_r, list(CASOS.keys()), active=0)
        for rot in self.radio.labels:
            rot.set_fontsize(10)
        self.radio.on_clicked(self._ao_escolher_caso)

        self.btn_pausa = Button(self.fig.add_axes([0.38, 0.14, 0.12, 0.055]),
                                "Pausar")
        self.btn_pausa.on_clicked(self._alternar_pausa)
        self.btn_reinicia = Button(self.fig.add_axes([0.38, 0.07, 0.12, 0.055]),
                                   "Reiniciar")
        self.btn_reinicia.on_clicked(self.reiniciar)
        self.fig.text(0.05, 0.025, f"Animação em câmera lenta ({ESCALA_TEMPO}×). "
                      "Mexer em qualquer controle reinicia a simulação.",
                      fontsize=8, color="#555555")

    def _ao_mudar(self, _valor):
        if not self._ajustando:
            self.reiniciar()

    def _ao_escolher_caso(self, rotulo):
        m1, m2, v1 = CASOS[rotulo]
        # _ajustando evita reiniciar a cada slider alterado por código
        self._ajustando = True
        self.s_m1.set_val(m1)
        self.s_m2.set_val(m2)
        self.s_v1.set_val(v1)
        self.s_v2.set_val(0.0)
        self._ajustando = False
        self.reiniciar()

    def _alternar_pausa(self, _evento=None):
        self.pausado = not self.pausado
        self.btn_pausa.label.set_text("Continuar" if self.pausado else "Pausar")
        self.fig.canvas.draw_idle()

    def reiniciar(self, _evento=None):
        self.cen = c = calcular_cenario(self.s_m1.val, self.s_m2.val,
                                        self.s_v1.val, self.s_v2.val,
                                        self.s_dt.val)
        self.t = 0.0
        self.pausado = False
        self.btn_pausa.label.set_text("Pausar")
        self.r1.set_width(c["w1"]); self.r1.set_height(c["w1"])
        self.r2.set_width(c["w2"]); self.r2.set_height(c["w2"])
        self.info1.set_text(f"Bloco 1  (m1 = {c['m1']:.1f} kg)")
        self.info2.set_text(f"Bloco 2  (m2 = {c['m2']:.1f} kg)")

        self._atualizar_barras(self.ax_p, self.bar_p, self.rot_p,
                               [c["p_antes"], c["p_depois"]])
        self._atualizar_barras(self.ax_k, self.bar_k, self.rot_k,
                               [c["K_antes"], c["K_depois"]])
        self._desenhar()

    @staticmethod
    def _atualizar_barras(ax, barras, rotulos, valores):
        ymin = min(0.0, min(valores)) * 1.3
        ymax = max(0.0, max(valores)) * 1.3
        if ymax - ymin < 1e-9:
            ymin, ymax = -1.0, 1.0
        ax.set_ylim(ymin, ymax)
        for b, r, v in zip(barras, rotulos, valores):
            b.set_height(v)
            r.set_position((b.get_x() + b.get_width() / 2, v))
            r.set_text(formata(v, sinal=False))
            r.set_va("bottom" if v >= 0 else "top")

    def atualizar(self, _quadro):
        c = self.cen
        if not self.pausado and self.t < c["t_fim"]:
            self.t = min(self.t + PASSO_T, c["t_fim"])
        self._desenhar()

    def _desenhar(self):
        c, t = self.cen, self.t
        x1, x2, va, vb = estado_em(c, t)
        self.r1.set_xy((x1 - c["w1"] / 2, 0))
        self.r2.set_xy((x2 - c["w2"] / 2, 0))
        self.n1.set_position((x1, c["w1"] / 2))
        self.n2.set_position((x2, c["w2"] / 2))

        depois = c["colide"] and t >= c["tc"]
        if not c["colide"]:
            estado = "sem colisão"
        elif depois:
            estado = "depois do choque"
        else:
            estado = "antes do choque"
        self.status.set_text(f"{estado}   (t = {t:.2f} s)")

        if not c["colide"]:
            self.aviso.set_text("Sem colisão: v1 ≤ v2 (os blocos não se aproximam)")
        elif c["fora"]:
            self.aviso.set_text(f"A colisão (t_c = {c['tc']:.1f} s) "
                                "ocorre fora da janela")
        elif t >= c["t_fim"]:
            self.aviso.set_text("Fim da simulação - clique em Reiniciar")
        else:
            self.aviso.set_text("")

        # barras "depois" esmaecidas até o choque
        alfa = 1.0 if (depois or not c["colide"]) else 0.35
        self.bar_p[1].set_alpha(alfa)
        self.bar_k[1].set_alpha(alfa)

        self.painel.set_text(self._texto_painel(t, va, vb, estado))

    def _texto_painel(self, t, va, vb, estado):
        c, f = self.cen, formata
        linha = "─" * 46
        L = ["RESULTADOS - colisão elástica (SI)", linha,
             f"{'':<15}{'antes':>14}{'depois':>14}",
             f"{'v1  (m/s)':<15}{f(c['v1']):>14}{f(c['u1']):>14}",
             f"{'v2  (m/s)':<15}{f(c['v2']):>14}{f(c['u2']):>14}",
             f"{'p   (kg·m/s)':<15}{f(c['p_antes']):>14}{f(c['p_depois']):>14}",
             f"{'K   (J)':<15}{f(c['K_antes'], sinal=False):>14}"
             f"{f(c['K_depois'], sinal=False):>14}",
             f"{'ΔK  (J)':<15}{'':>14}{f(c['dK']):>14}", linha]

        if c["colide"]:
            L += [f"Instante do choque: t_c = {c['tc']:.3f} s",
                  f"Impulso e força média (Δt = {c['dt_ms']:.0f} ms)",
                  f"J1 = {f(c['J1'])} N·s    J2 = {f(c['J2'])} N·s",
                  f"J1 + J2 = {f(c['J1'] + c['J2'])} N·s",
                  f"F_med = {c['F_med']:.2f} N",
                  f"u2 - u1 = {f(c['u2'] - c['u1'])}   "
                  f"-(v2 - v1) = {f(-(c['v2'] - c['v1']))}"]
        else:
            L += ["SEM COLISÃO (v1 <= v2): os blocos não se",
                  "encontram; velocidades inalteradas,",
                  "J1 = J2 = 0 e F_med = 0."]

        pa = momento(c["m1"], va, c["m2"], vb)
        ka = cinetica(c["m1"], va, c["m2"], vb)
        L += [linha, f"AGORA: t = {t:.2f} s  ({estado})",
              f"v1 = {f(va)} m/s    v2 = {f(vb)} m/s",
              f"p = {f(pa)} kg·m/s    K = {f(ka, sinal=False)} J"]
        return "\n".join(L)


def main():
    interface = Interface()
    plt.show()
    return interface


if __name__ == "__main__":
    main()
