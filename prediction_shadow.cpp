//
// このコードは Hsin-Yuan Huang (https://momohuang.github.io/)
// によって作成されました。 詳細は以下の論文を参照してください:
//  "Predicting Many Properties of a Quantum System from Very Few Measurements"
//  (非常に少ない測定から量子系の多くの特性を予測する)
//
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <utility>
#include <vector>

using namespace std;

int system_size = -1;
int number_of_observables;

double renyi_sum_of_binary_outcome[100000000];
double renyi_number_of_outcomes[100000000];

//
// 以下の関数はファイル: observable_file_name を読み込み、
// [observables] と [observables_acting_on_ith_qubit] を更新します。
//
vector<vector<pair<int, int>>> observables; // 予測する観測量 (Pauli observable)
vector<vector<vector<int>>> observables_acting_on_ith_qubit;
void read_all_observables(char *observable_file_name) {
  ifstream observable_fstream;
  observable_fstream.open(observable_file_name, ifstream::in);

  if (observable_fstream.fail()) {
    fprintf(stderr,
            "\n====\nError: 入力ファイル \"%s\" が存在しません。\n====\n",
            observable_file_name);
    exit(-1);
  }

  // システムサイズを読み込む
  int system_size_observable;
  observable_fstream >> system_size_observable;
  if (system_size == -1)
    system_size = system_size_observable;

  // 以下の数値を初期化 (後で変更される)
  number_of_observables = 0;

  //
  // observables_acting_on_ith_qubit の構造:
  //   observables_acting_on_ith_qubit[ith_qubit][0 or 1 or 2]
  //     ith_qubit に対して X (0), Y (1), Z (2)
  //     を適用する観測量のインデックスリストを返します
  //
  observables_acting_on_ith_qubit.clear();
  vector<int> single_list;
  vector<vector<int>> pauli_list;
  pauli_list.push_back(single_list);
  pauli_list.push_back(single_list);
  pauli_list.push_back(single_list);
  for (int i = 0; i < system_size; i++) {
    observables_acting_on_ith_qubit.push_back(pauli_list);
  }

  // 局所観測量を行ごとに読み込む
  string line;
  int observable_counter = 0;
  while (getline(observable_fstream, line)) {
    if (line == "\n" || line == "")
      continue;
    istringstream single_line_stream(line);

    int k_local;
    single_line_stream >> k_local;

    vector<pair<int, int>> ith_observable;

    for (int k = 0; k < k_local; k++) {
      char pauli_observable[5];
      int position_of_pauli;
      single_line_stream >> pauli_observable >> position_of_pauli;

      assert(pauli_observable[0] == 'X' || pauli_observable[0] == 'Y' ||
             pauli_observable[0] == 'Z');

      int pauli_encoding = pauli_observable[0] - 'X'; // X -> 0, Y -> 1, Z -> 2

      observables_acting_on_ith_qubit[position_of_pauli][pauli_encoding]
          .push_back(observable_counter);
      ith_observable.push_back(make_pair(position_of_pauli, pauli_encoding));
    }

    observables.push_back(ith_observable);
    observable_counter++;
  }
  number_of_observables = observable_counter;
  observable_fstream.close();

  return;
}

//
// 以下の関数はファイル: subsystem_file_name を読み込み、
// [subsystems] を更新します。
//
vector<vector<int>> subsystems; // エントロピーを予測する部分系
void read_all_subsystems(char *subsystem_file_name) {
  ifstream subsystem_fstream;
  subsystem_fstream.open(subsystem_file_name, ifstream::in);

  if (subsystem_fstream.fail()) {
    fprintf(stderr,
            "\n====\nError: 入力ファイル \"%s\" が存在しません。\n====\n",
            subsystem_file_name);
    exit(-1);
  }

  // システムサイズを読み込む
  int system_size_subsystem;
  subsystem_fstream >> system_size_subsystem;
  if (system_size == -1)
    system_size = system_size_subsystem;

  // 局所観測量を行ごとに読み込む
  string line;
  int observable_counter = 0;
  while (getline(subsystem_fstream, line)) {
    if (line == "\n" || line == "")
      continue;
    istringstream single_line_stream(line);

    int k_local;
    single_line_stream >> k_local;

    vector<int> ith_subsystem;

    for (int k = 0; k < k_local; k++) {
      int position_of_the_qubit;
      single_line_stream >> position_of_the_qubit;
      ith_subsystem.push_back(position_of_the_qubit);
    }

    subsystems.push_back(ith_subsystem);
  }
  subsystem_fstream.close();

  return;
}

//
// 以下の関数はファイル: measurement_file_name を読み込み、
// [observables] と [observables_acting_on_ith_qubit] を更新します。
//
vector<vector<int>> measurement_pauli_basis;
vector<vector<int>> measurement_binary_outcome;
void read_all_measurements(char *measurement_file_name) {
  ifstream measurement_fstream;
  measurement_fstream.open(measurement_file_name, ifstream::in);

  if (measurement_fstream.fail()) {
    fprintf(stderr,
            "\n====\nError: 入力ファイル \"%s\" が存在しません。\n====\n",
            measurement_file_name);
    exit(-1);
  }

  // システムサイズを読み込む
  int system_size_measurement;
  measurement_fstream >> system_size_measurement;
  if (system_size == -1)
    system_size = system_size_measurement;
  if (system_size_measurement != system_size) {
    fprintf(stderr, "\n====\nError: システムサイズが一致しません。\n====\n");
    exit(-1);
  }

  // 測定結果を行ごとに読み込む
  string line;
  int measurement_counter = 0;
  while (getline(measurement_fstream, line)) {
    if (line == "\n" || line == "")
      continue;
    istringstream single_line_stream(line);

    vector<int> empty_list;
    measurement_pauli_basis.push_back(empty_list);
    measurement_binary_outcome.push_back(empty_list);
    for (int ith_qubit = 0; ith_qubit < system_size; ith_qubit++) {
      char pauli[10];
      int binary_outcome;
      single_line_stream >> pauli >> binary_outcome;
      assert(binary_outcome == 1 || binary_outcome == -1);

      measurement_pauli_basis[measurement_counter].push_back(pauli[0] - 'X');
      measurement_binary_outcome[measurement_counter].push_back(binary_outcome);
    }

    measurement_counter++;
  }
}

//
// 以下の関数はこのプログラムの使用法を表示します。
//
void print_usage() {
  fprintf(stderr, "使用法:\n");
  fprintf(stderr,
          "./prediction_shadow -o [measurement.txt] [observable.txt]\n");
  fprintf(stderr, "    このオプションは局所観測量の期待値を予測します。\n");
  fprintf(stderr, "    [observable.txt] "
                  "で指定された各局所観測量について、予測値を出力します。\n");
  fprintf(stderr, "<または>\n");
  fprintf(stderr, "./prediction_shadow -e [subsystem.txt] [measurement.txt]\n");
  fprintf(stderr, "    このオプションは Renyi (レニー) "
                  "エンタングルメントエントロピーを予測します。\n");
  fprintf(
      stderr,
      "    [subsystem.txt] "
      "で指定された各部分系について、予測されたエントロピーを出力します。\n");
  return;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    print_usage();
    return -1;
  }

  //
  // 局所観測量の予測を実行
  //
  if (strcmp(argv[1], "-o") == 0) {
    read_all_measurements(argv[2]);
    read_all_observables(argv[3]);

    // すべての観測量について、
    // 現在の測定繰り返しにおいて、その観測量を測定するために
    // いくつのパウリ演算子が一致する必要があるか
    vector<int> how_many_pauli_to_match;
    how_many_pauli_to_match.resize(number_of_observables);

    // すべての観測量について、
    // この単一量子ビット測定までの累積結果を保存
    vector<int> cumulative_measurement;
    cumulative_measurement.resize(number_of_observables);

    // すべての観測量について、
    // それが何回測定されたかを保存 (マッチした場合のみカウント)
    vector<int> number_of_measurements;
    number_of_measurements.resize(number_of_observables);

    // すべての観測量について、
    // 測定結果の合計を保存
    vector<int> sum_of_measurement_results;
    sum_of_measurement_results.resize(number_of_observables);

    // 測定データを走査:
    //    measurement_pauli_basis, measurement_binary_outcome
    // 局所観測量を計算するために
    for (int t = 0; t < (int)measurement_pauli_basis.size(); t++) {
      for (int i = 0; i < (int)observables.size(); i++) {
        how_many_pauli_to_match[i] =
            observables[i].size();     // k-local 観測量の場合は k で初期化
        cumulative_measurement[i] = 1; // 1 で初期化
      }

      for (int ith_qubit = 0; ith_qubit < system_size; ith_qubit++) {
        int pauli = measurement_pauli_basis[t][ith_qubit];
        int binary_outcome = measurement_binary_outcome[t][ith_qubit];
        for (int i : observables_acting_on_ith_qubit[ith_qubit][pauli]) {
          how_many_pauli_to_match[i]--;
          cumulative_measurement[i] *= binary_outcome;
        }
      }

      for (int i = 0; i < (int)observables.size(); i++) {
        if (how_many_pauli_to_match[i] == 0) {
          number_of_measurements[i]++;
          sum_of_measurement_results[i] += cumulative_measurement[i];
        }
      }
    }

    for (int i = 0; i < (int)observables.size(); i++) {
      if (number_of_measurements[i] == 0) {
        fprintf(stderr, "%d-th Observable is not measured at all\n", i + 1);
        printf("0\n");
      } else
        printf("%f\n",
               1.0 * sum_of_measurement_results[i] / number_of_measurements[i]);
    }
  }
  //
  // エンタングルメントエントロピーの予測を実行
  //
  else if (strcmp(argv[1], "-e") == 0) {
    read_all_measurements(argv[2]);
    read_all_subsystems(argv[3]);

    for (int s = 0; s < (int)subsystems.size(); s++) {
      int subsystem_size = (int)subsystems[s].size();

      for (int c = 0; c < (1 << (2 * subsystem_size)); c++) {
        renyi_sum_of_binary_outcome[c] = 0;
        renyi_number_of_outcomes[c] = 0;
      }

      for (int t = 0; t < (int)measurement_pauli_basis.size(); t++) {
        long long encoding = 0, cumulative_outcome = 1;

        renyi_sum_of_binary_outcome[0] += 1;
        renyi_number_of_outcomes[0] += 1;

        // グレイコード (Gray code) を使用して、すべての 2^n
        // 個の可能な結果を反復処理
        for (long long b = 1; b < (1 << subsystem_size); b++) {
          long long change_i = __builtin_ctzll(b);
          long long index_in_original_system = subsystems[s][change_i];

          cumulative_outcome *=
              measurement_binary_outcome[t][index_in_original_system];
          encoding ^= (measurement_pauli_basis[t][index_in_original_system] + 1)
                      << (2LL * change_i);

          renyi_sum_of_binary_outcome[encoding] += cumulative_outcome;
          renyi_number_of_outcomes[encoding] += 1;
        }
      }
      vector<int> level_cnt(2 * subsystem_size, 0);
      vector<int> level_ttl(2 * subsystem_size, 0);

      for (long long c = 0; c < (1 << (2 * subsystem_size)); c++) {
        int nonId = 0;
        for (int i = 0; i < subsystem_size; i++) {
          nonId += ((c >> (2 * i)) & 3) != 0;
        }
        if (renyi_number_of_outcomes[c] >= 2)
          level_cnt[nonId]++;
        level_ttl[nonId]++;
      }

      double predicted_entropy = 0;
      for (long long c = 0; c < (1 << (2 * subsystem_size)); c++) {
        if (renyi_number_of_outcomes[c] <= 1)
          continue;

        int nonId = 0;
        for (int i = 0; i < subsystem_size; i++)
          nonId += ((c >> (2 * i)) & 3) != 0;

        predicted_entropy +=
            ((double)1.0) /
            (renyi_number_of_outcomes[c] * (renyi_number_of_outcomes[c] - 1)) *
            (renyi_sum_of_binary_outcome[c] * renyi_sum_of_binary_outcome[c] -
             renyi_number_of_outcomes[c]) /
            (1LL << subsystem_size) * level_ttl[nonId] / level_cnt[nonId];
      }

      printf("%f\n", -1.0 * log2(min(max(predicted_entropy,
                                         1.0 / pow(2.0, subsystem_size)),
                                     1.0 - 1e-9)));
    }
  }
  //
  // 上記のいずれにも該当しない (入力が無効)
  //
  else {
    print_usage();
    return -1;
  }
}
