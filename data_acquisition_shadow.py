#
# このコードは Hsin-Yuan Huang (https://momohuang.github.io/) によって作成されました。
# 詳細は以下の論文を参照してください:
#  "Predicting Many Properties of a Quantum System from Very Few Measurements"
#  (非常に少ない測定から量子系の多くの特性を予測する)
# この Python バージョンは C++ バージョンよりも低速です。(コードの最適化が少ないため)
# しかし、理解しやすく、拡張もしやすいはずです。
#
import sys
import random
import math

def randomized_classical_shadow(num_total_measurements, system_size):
    #
    # 古典シャドウのランダム化バージョンの実装
    #
    #    num_total_measurements: int, 測定ラウンドの総数 (T)
    #    system_size: int, 量子系の量子ビット数 (n)
    #
    measurement_procedure = []
    for t in range(num_total_measurements):
        # 各量子ビットに対してランダムにパウリ基底 (X, Y, Z) を選択
        single_round_measurement = [random.choice(["X", "Y", "Z"]) for i in range(system_size)]
        measurement_procedure.append(single_round_measurement)
    return measurement_procedure

def derandomized_classical_shadow(all_observables, num_of_measurements_per_observable, system_size, weight=None):
    #
    # 古典シャドウの非ランダム化バージョンの実装
    #
    #     all_observables: パウリ観測量のリスト。各パウリ観測量はタプルのリストです。
    #                       ("X", position) または ("Y", position) または ("Z", position) の形式
    #     num_of_measurements_per_observable: int, 各観測量に対する測定回数
    #     system_size: int, 量子系の量子ビット数
    #     weight: None または 各観測量に対する係数のリスト
    #             None -- このパラメータを無視
    #             リスト -- 対応する重みで各観測量の測定回数を変更
    #
    if weight is None:
        weight = [1.0] * len(all_observables)
    assert(len(weight) == len(all_observables))

    sum_log_value = 0
    sum_cnt = 0

    def cost_function(num_of_measurements_so_far, num_of_matches_needed_in_this_round, shift = 0):
        eta = 0.9 # 変更可能なハイパーパラメータ
        nu = 1 - math.exp(-eta / 2)

        nonlocal sum_log_value
        nonlocal sum_cnt

        cost = 0
        for i, zipitem in enumerate(zip(num_of_measurements_so_far, num_of_matches_needed_in_this_round)):
            measurement_so_far, matches_needed = zipitem
            # 既に十分な回数測定されている場合はスキップ
            if num_of_measurements_so_far[i] >= math.floor(weight[i] * num_of_measurements_per_observable):
                continue

            if system_size < matches_needed:
                V = eta / 2 * measurement_so_far
            else:
                V = eta / 2 * measurement_so_far - math.log(1 - nu / (3 ** matches_needed))
            cost += math.exp(-V / weight[i] - shift)

            sum_log_value += V / weight[i]
            sum_cnt += 1

        return cost

    def match_up(qubit_i, dice_roll_pauli, single_observable):
        for pauli, pos in single_observable:
            if pos != qubit_i:
                continue
            else:
                if pauli != dice_roll_pauli:
                    return -1 # 不一致
                else:
                    return 1 # 一致
        return 0 # 無関係

    num_of_measurements_so_far = [0] * len(all_observables)
    measurement_procedure = []

    for repetition in range(num_of_measurements_per_observable * len(all_observables)):
        # "system_size" 個の量子ビットに対する並列測定の1ラウンド
        num_of_matches_needed_in_this_round = [len(P) for P in all_observables]
        single_round_measurement = []

        shift = sum_log_value / sum_cnt if sum_cnt > 0 else 0;
        sum_log_value = 0.0
        sum_cnt = 0

        for qubit_i in range(system_size):
            cost_of_outcomes = dict([("X", 0), ("Y", 0), ("Z", 0)])

            for dice_roll_pauli in ["X", "Y", "Z"]:
                # 測定結果が "dice_roll_pauli" になると仮定
                for i, single_observable in enumerate(all_observables):
                    result = match_up(qubit_i, dice_roll_pauli, single_observable)
                    if result == -1:
                        num_of_matches_needed_in_this_round[i] += 100 * (system_size+10) # 測定不可能
                    if result == 1:
                        num_of_matches_needed_in_this_round[i] -= 1 # Pauli X/Y/Z が1つ一致

                cost_of_outcomes[dice_roll_pauli] = cost_function(num_of_measurements_so_far, num_of_matches_needed_in_this_round, shift=shift)

                # 仮定を元に戻す
                for i, single_observable in enumerate(all_observables):
                    result = match_up(qubit_i, dice_roll_pauli, single_observable)
                    if result == -1:
                        num_of_matches_needed_in_this_round[i] -= 100 * (system_size+10) # 測定不可能
                    if result == 1:
                        num_of_matches_needed_in_this_round[i] += 1 # Pauli X/Y/Z が1つ一致

            for dice_roll_pauli in ["X", "Y", "Z"]:
                if min(cost_of_outcomes.values()) < cost_of_outcomes[dice_roll_pauli]:
                    continue
                # コストが最小の測定基底をここで選択
                single_round_measurement.append(dice_roll_pauli)
                for i, single_observable in enumerate(all_observables):
                    result = match_up(qubit_i, dice_roll_pauli, single_observable)
                    if result == -1:
                        num_of_matches_needed_in_this_round[i] += 100 * (system_size+10) # 測定不可能
                    if result == 1:
                        num_of_matches_needed_in_this_round[i] -= 1 # Pauli X/Y/Z が1つ一致
                break

        measurement_procedure.append(single_round_measurement)

        for i, single_observable in enumerate(all_observables):
            if num_of_matches_needed_in_this_round[i] == 0: # すべての量子ビットの測定が完了
                num_of_measurements_so_far[i] += 1

        success = 0
        for i, single_observable in enumerate(all_observables):
            if num_of_measurements_so_far[i] >= math.floor(weight[i] * num_of_measurements_per_observable):
                success += 1

        if success == len(all_observables):
            break

    return measurement_procedure

#
# 以下のコードは、コマンドラインインターフェースを介してこのコードを実行する場合にのみ使用されます
#
if __name__ == "__main__":
    def print_usage():
        print("Usage:\n", file=sys.stderr)
        print("python3 data_acquisition_shadow -d [観測量ごとの測定回数] [observable.txt]", file=sys.stderr)
        print("    これは古典シャドウの非ランダム化 (derandomized) バージョンです。", file=sys.stderr)
        print("    [observable.txt] にあるすべての観測量を少なくとも [観測量ごとの測定回数] 回測定するための", file=sys.stderr)
        print("    パウリ測定のリストを出力します。", file=sys.stderr)
        print("<or>\n", file=sys.stderr)
        print("python3 data_acquisition_shadow -r [総測定回数] [システムサイズ]", file=sys.stderr)
        print("    これは古典シャドウのランダム化 (randomized) バージョンです。", file=sys.stderr)
        print("    指定された [システムサイズ] に対して、合計 [総測定回数] 回の繰り返しのための", file=sys.stderr)
        print("    パウリ測定のリストを出力します。", file=sys.stderr)

    if len(sys.argv) != 4:
        print_usage()
    if sys.argv[1] == "-d":
        with open(sys.argv[3]) as f:
            content = f.readlines()
        system_size = int(content[0])

        all_observables = []
        for line in content[1:]:
            one_observable = []
            for pauli_XYZ, position in zip(line.split(" ")[1::2], line.split(" ")[2::2]):
                one_observable.append((pauli_XYZ, int(position)))
            all_observables.append(one_observable)
        measurement_procedure = derandomized_classical_shadow(all_observables, int(sys.argv[2]), system_size)
    elif sys.argv[1] == "-r":
        measurement_procedure = randomized_classical_shadow(int(sys.argv[2]), int(sys.argv[3]))
    else:
        print_usage()

    for single_round_measurement in measurement_procedure:
        print(" ".join(single_round_measurement))
