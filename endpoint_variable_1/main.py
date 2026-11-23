import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))  # 공통 모듈
import numpy as np
import matplotlib.pyplot as plt

import nodes
import draw_basis as draw
from optimal_nodes_maximum import optimal_node


N_MIN = 4
N_MAX = 35

if __name__ == "__main__":

    # GL, GLo, ESC는 생성 함수를 넘겨 선택한 n에 맞춰 기저를 새로 계산한다.
    # newton은 optimal_nodes_maximum.py에 나열된 n = 4 ~ 35의 Gauss-Newton 결과 노드를 쓴다.
    basis_sources = {
        "GL": nodes.GL_nodes,
        "GLo": nodes.GLo_nodes,
        "ESC": nodes.ESC_nodes,
        "newton": {int(name[1:]): np.array(node) for name, node in optimal_node.items()},
    }

    ##################################################################################################################
    # 1️⃣ n = 4 ~ 35에서 GL / GLo / ESC / Gauss-Newton 결과 노드의 기저 비교
    # 슬라이더(또는 좌우 방향키)로 n을, 라디오 버튼으로 노드 종류를 고른다.
    fig, ax, basis_selector, degree_slider = draw.draw_basis_selector(
        basis_sources,
        n_initial=N_MIN,
        n_min=N_MIN,
        n_max=N_MAX,
        initial_basis="newton",
    )

    ##################################################################################################################
    # 2️⃣ 등간격 / GLo / ESC / Gauss-Newton 최적 노드의 [-1, 1] 위 배치
    # 체크 박스로 비교할 노드를 여러 개 고르고, 슬라이더(또는 좌우 방향키)로 n을 고른다.
    # 최적 노드는 optimal_nodes_maximum.py에 있는 n = 4 ~ 55까지 있다.
    layout_sources = {
        "uniform": nodes.uniform_nodes,
        "GLo": nodes.GLo_nodes,
        "ESC": nodes.ESC_nodes,
        "newton": {int(name[1:]): node for name, node in optimal_node.items()},
    }
    layout_fig, layout_ax, layout_checks, layout_slider = draw.draw_node_layout(
        layout_sources,
        n_initial=N_MIN,
        n_min=N_MIN,
        n_max=max(layout_sources["newton"]),
    )

    plt.show()
