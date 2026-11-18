# ADO-3 - Colisões, Impulso e Conservação do Momento Linear

Física: Práticas Extensivas - Prof. Ettore Baldini-Neto

Simulação interativa de uma colisão frontal **elástica** entre dois blocos (Opção A),
com animação, gráficos de momento e energia e painel de resultados.

## Como executar
```bash
pip install matplotlib numpy
python ADO3_FISICA.py
```

## O que o programa mostra
- Animação dos dois blocos em MRU antes e depois do choque
- Sliders para m1, m2, v1, v2 e Δt de contato
- Menu com os 3 casos (m1 = m2, m1 < m2, m1 > m2), botões Pausar e Reiniciar
- Painel com velocidades finais, p, K, ΔK, impulsos J1 e J2 e força média
- Gráficos de barras de p e K, antes e depois
- Aviso na tela quando não há colisão (v1 ≤ v2)
