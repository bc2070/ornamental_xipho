import pandas as pd
import os
import sys


def parse_state(state):
    """将 '2,0,0' 转换为 (2,0,0)。"""
    return tuple(int(x) for x in state.split(','))


def detect_ancestry_count(states):
    """根据 genotype state 判断 ancestry 数量。"""
    if not states:
        return 0

    return len(parse_state(states[0]))


def calculate_contribution(df, states, ancestry_count):
    """
    根据整个 posterior 文件计算每个 ancestry 的贡献比例。

    对三祖先：

        A = 2*S200 + S110 + S101
        B = S110 + 2*S020 + S011
        C = S101 + S011 + 2*S002

    最后除以所有 ancestry allele contribution 总和。
    """

    contribution = [0.0] * ancestry_count

    for state in states:

        genotype = parse_state(state)

        posterior_sum = df[state].sum()

        for i in range(ancestry_count):

            contribution[i] += (
                genotype[i] * posterior_sum
            )

    total = sum(contribution)

    if total == 0:

        return [0.0] * ancestry_count

    return [
        x / total
        for x in contribution
    ]


def make_state_mapping(states, ranking, ancestry_count):
    """
    将三祖先 genotype 转换成标准双祖先 genotype。

    ranking:
        例如 [2, 0, 1]

    表示：

        Rank1 = 原始 ancestry 3
        Rank2 = 原始 ancestry 1
        Rank3 = 原始 ancestry 2

    只保留 Rank1 和 Rank2。

    Rank3 相关 genotype -> NA
    """

    top1 = ranking[0]
    top2 = ranking[1]

    mapping = {}

    for state in states:

        genotype = parse_state(state)

        # --------------------------------------------------
        # 检查是否包含 Rank3
        # --------------------------------------------------

        third_count = sum(
            genotype[i]
            for i in range(ancestry_count)
            if i not in [top1, top2]
        )

        if third_count > 0:

            mapping[state] = "NA"

            continue

        # --------------------------------------------------
        # Top1 / Top2 genotype
        # --------------------------------------------------

        count1 = genotype[top1]
        count2 = genotype[top2]

        if count1 == 2 and count2 == 0:

            mapping[state] = "2,0"

        elif count1 == 1 and count2 == 1:

            mapping[state] = "1,1"

        elif count1 == 0 and count2 == 2:

            mapping[state] = "0,2"

        else:

            mapping[state] = "NA"

    return mapping


def transform_posterior(
    input_file,
    output_file,
    ranking_file
):

    # ======================================================
    # Read posterior
    # ======================================================

    df = pd.read_csv(
        input_file,
        sep=r'\s+'
    )

    if len(df.columns) < 3:

        raise ValueError(
            f"{input_file}: posterior 格式错误"
        )

    metadata_columns = [
        df.columns[0],
        df.columns[1]
    ]

    states = list(
        df.columns[2:]
    )

    ancestry_count = detect_ancestry_count(
        states
    )

    print()
    print("=" * 60)
    print(
        f"Processing: {os.path.basename(input_file)}"
    )
    print("=" * 60)

    print(
        f"Detected ancestry count: {ancestry_count}"
    )

    # ======================================================
    # Two ancestry
    # ======================================================

    if ancestry_count == 2:

        print(
            "Status: TWO_ANCESTRY"
        )

        df.to_csv(
            output_file,
            sep='\t',
            index=False
        )

        with open(
            ranking_file,
            'w'
        ) as f:

            f.write(
                "ANCESTRY_COUNT\t2\n"
            )

            f.write(
                "STATUS\tTWO_ANCESTRY\n"
            )

            f.write(
                "NO_TRANSFORMATION\tTRUE\n"
            )

        print(
            f"[OK] Output: {output_file}"
        )

        return

    # ======================================================
    # Only support 2 or 3 ancestry
    # ======================================================

    if ancestry_count != 3:

        raise ValueError(
            f"{input_file}: "
            f"目前只支持 2 或 3 ancestry，"
            f"检测到 {ancestry_count}"
        )

    # ======================================================
    # Calculate contribution
    # ======================================================

    contribution = calculate_contribution(
        df,
        states,
        ancestry_count
    )

    # ======================================================
    # Rank ancestry
    #
    # 注意：
    # 这里按照 contribution 排序，
    # 不是按照 posterior column 顺序。
    # ======================================================

    ranking = sorted(
        range(ancestry_count),
        key=lambda i: contribution[i],
        reverse=True
    )

    print()
    print("Ancestry contribution:")

    for i in range(ancestry_count):

        print(
            f"  Original ancestry {i + 1}: "
            f"{contribution[i]:.8f}"
        )

    print()
    print("Ancestry ranking:")

    for rank, idx in enumerate(
        ranking,
        start=1
    ):

        print(
            f"  Rank{rank} = "
            f"Original ancestry {idx + 1} "
            f"({contribution[idx]:.8f})"
        )

    # ======================================================
    # State mapping
    # ======================================================

    state_mapping = make_state_mapping(
        states,
        ranking,
        ancestry_count
    )

    print()
    print("State mapping:")

    for old_state in states:

        print(
            f"  {old_state} -> "
            f"{state_mapping[old_state]}"
        )

    # ======================================================
    # Create standard two-ancestry posterior
    # ======================================================

    output = df[
        metadata_columns
    ].copy()

    output["2,0"] = 0.0
    output["1,1"] = 0.0
    output["0,2"] = 0.0

    # ======================================================
    # Project three ancestry posterior
    #
    # Rank3-related genotype is ignored.
    # ======================================================

    for old_state in states:

        new_state = state_mapping[
            old_state
        ]

        if new_state in [
            "2,0",
            "1,1",
            "0,2"
        ]:

            output[new_state] += (
                df[old_state]
            )

    # ======================================================
    # Renormalize Top1 / Top2 posterior
    #
    # 例如：
    #
    # 2,0 = 0.80
    # 1,1 = 0.10
    # 0,2 = 0.03
    #
    # total = 0.93
    #
    # 最终：
    #
    # 2,0 = 0.80 / 0.93
    # 1,1 = 0.10 / 0.93
    # 0,2 = 0.03 / 0.93
    #
    # 使三列重新满足：
    #
    # P(2,0) + P(1,1) + P(0,2) = 1
    #
    # 对存在 Top1/Top2 posterior 的位置成立。
    # ======================================================

    two_ancestry_states = [
        "2,0",
        "1,1",
        "0,2"
    ]

    two_ancestry_total = (
        output["2,0"]
        + output["1,1"]
        + output["0,2"]
    )

    # 避免除以 0
    valid = two_ancestry_total > 0

    for state in two_ancestry_states:

        output.loc[valid, state] = (
            output.loc[valid, state]
            / two_ancestry_total[valid]
        )

    # 如果三个概率都为 0，则保持为 0
    output.loc[
        ~valid,
        two_ancestry_states
    ] = 0.0

    # ======================================================
    # Save transformed posterior
    # ======================================================

    output.to_csv(
        output_file,
        sep='\t',
        index=False
    )

    # ======================================================
    # Save ancestry ranking / mapping information
    # ======================================================

    with open(
        ranking_file,
        'w'
    ) as f:

        f.write(
            "ANCESTRY_COUNT\t3\n"
        )

        f.write(
            "STATUS\tTHREE_TO_TWO\n"
        )

        f.write(
            "TOP_ANCESTRIES\t2\n"
        )

        f.write(
            "MINOR_ANCESTRY_IGNORED\tTRUE\n"
        )

        f.write(
            "\nANCESTRY_RANKING\n"
        )

        for rank, idx in enumerate(
            ranking,
            start=1
        ):

            f.write(
                f"RANK{rank}\t"
                f"ORIGINAL_INDEX={idx + 1}\t"
                f"CONTRIBUTION="
                f"{contribution[idx]:.12f}\n"
            )

        f.write(
            "\nSTATE_MAPPING\n"
        )

        for old_state in states:

            f.write(
                f"{old_state}\t"
                f"{state_mapping[old_state]}\n"
            )

    # ======================================================
    # Summary
    # ======================================================

    print()
    print(
        f"[OK] Output posterior: "
        f"{output_file}"
    )

    print(
        f"[OK] Ranking file: "
        f"{ranking_file}"
    )

    print(
        "Final posterior columns:"
    )

    print(
        "  chrom position 2,0 1,1 0,2"
    )

    print()


def main():

    if len(sys.argv) != 4:

        print(
            "Usage:"
        )

        print(
            "python3 00.prepare_posterior.py "
            "input.posterior "
            "output.posterior "
            "ranking.txt"
        )

        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ranking_file = sys.argv[3]

    transform_posterior(
        input_file,
        output_file,
        ranking_file
    )


if __name__ == "__main__":

    main()
