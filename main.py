import numpy as np
import matplotlib.pyplot as plt

import generate_nodes as nodes
import draw_basis as draw
from optimal_nodes_maximum import optimal_node


N_MIN = 4
N_MAX = 15


if __name__ == "__main__":

    # GL, GLo, ESC는 생성 함수를 넘겨 선택한 n에 맞춰 기저를 새로 계산한다.
    # compass는 optimal_nodes_maximum.py에 나열된 n = 4 ~ 15의 compass search 결과 노드를 쓴다.
    basis_sources = {
        "GL": nodes.GL_nodes,
        "GLo": nodes.GLo_nodes,
        "ESC": nodes.ESC_nodes,
        "compass": {int(name[1:]): np.array(node) for name, node in optimal_node.items()},
    }

    ##################################################################################################################
    # 1️⃣ n = 4 ~ 15에서 GL / GLo / ESC / compass search 결과 노드의 기저 비교
    # 슬라이더(또는 좌우 방향키)로 n을, 라디오 버튼으로 노드 종류를 고른다.
    fig, ax, basis_selector, degree_slider = draw.draw_basis_selector(
        basis_sources,
        n_initial=N_MIN,
        n_min=N_MIN,
        n_max=N_MAX,
        initial_basis="compass",
    )

    plt.show()
