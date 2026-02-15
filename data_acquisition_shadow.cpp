//
// このコードは Hsin-Yuan Huang (https://momohuang.github.io/)
// によって作成されました。 詳細は以下の論文を参照してください:
//  "Predicting Many Properties of a Quantum System from Very Few Measurements"
//  (非常に少ない測定から量子系の多くの特性を予測する)
//
#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <utility>
#include <vector>

using namespace std;
const int INF = 999999999; // 無限大として扱う非常に大きな数

double eta = 0.9; // チューニングが必要なハイパーパラメータ

int system_size;
int number_of_observables;
int number_of_measurements_per_observable;
int max_k_local;

//
// 以下の関数はファイル: observable_file_name を読み込み、
// [observables] と [observables_acting_on_ith_qubit] を更新します。
//
vector<vector<pair<int, int>>> observables; // 予測する観測量 (Pauli observable)
vector<vector<vector<int>>> observables_acting_on_ith_qubit;
vector<double> observables_weight;

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
  observable_fstream >> system_size;

  // 以下の数値を初期化 (後で変更される)
  max_k_local = 0;
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
    if (line == "\n")
      continue;
    istringstream single_line_stream(line);

    int k_local;
    single_line_stream >> k_local;
    max_k_local = max(max_k_local, k_local);

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

    double weight;
    int X = single_line_stream.rdbuf()->in_avail();
    if (X == 0)
      weight = 1.0;
    else
      single_line_stream >> weight;

    observables_weight.push_back(weight);

    observables.push_back(ith_observable);
    observable_counter++;
  }
  number_of_observables = observable_counter;
  observable_fstream.close();

  return;
}

//
// 以下の関数はこのプログラムの使用法を表示します。
//
void print_usage() {
  fprintf(stderr, "使用法:\n");
  fprintf(
      stderr,
      "./shadow_data_acquisition -d [観測量ごとの測定回数] [observable.txt]\n");
  fprintf(
      stderr,
      "    これは古典シャドウの非ランダム化 (derandomized) バージョンです。\n");
  fprintf(stderr, "    [observable.txt] にあるすべての観測量を少なくとも "
                  "[観測量ごとの測定回数] 回測定するための\n");
  fprintf(stderr, "    パウリ測定のリストを出力します。\n");
  fprintf(stderr, "<または>\n");
  fprintf(stderr,
          "./shadow_data_acquisition -r [総測定回数] [システムサイズ]\n");
  fprintf(stderr,
          "    これは古典シャドウのランダム化 (randomized) バージョンです。\n");
  fprintf(stderr, "    指定された [システムサイズ] に対して、合計 [総測定回数] "
                  "回の繰り返しのための\n");
  fprintf(stderr, "    パウリ測定のリストを出力します。\n");
  return;
}

//
// 以下の関数は乗法重み更新 (multiplicative weight update) 法を実行します。
// これは古典シャドウのランダムなパウリ測定を非ランダム化するために使用されます。
//
vector<double>
    log1ppow1o3k; // log1ppow1o3k[k] = log(1 + (e^(-eta / 2) - 1) / 3^k)
double sum_log_value = 0.0;
int sum_cnt = 0.0;
double
fail_prob_pessimistic(int cur_num_of_measurements, int how_many_pauli_to_match,
                      double weight,
                      double shift) { // "悲観的推定による失敗確率" を意味します
  double log1pp0 =
      (how_many_pauli_to_match < INF ? log1ppow1o3k[how_many_pauli_to_match]
                                     : 0.0);

  if (floor(weight * number_of_measurements_per_observable) <=
      cur_num_of_measurements)
    return 0;

  double log_value = -eta / 2 * cur_num_of_measurements + log1pp0;
  sum_log_value += (log_value / weight);
  sum_cnt++;
  return 2 * exp((log_value / weight) - shift);
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    print_usage();
    return -1;
  }

  //
  // 古典シャドウのランダム化バージョンを実行
  //
  if (strcmp(argv[1], "-r") == 0) {
    //
    // この実行のためのランダムシードを設定
    //
    srand(time(NULL));

    //
    // パラメータを読み込む
    //
    system_size = stoi(argv[3]);
    int number_of_total_measurements = stoi(argv[2]);
    char Pauli[] = {'X', 'Y', 'Z'};

    //
    // 古典シャドウのランダム化バージョン
    // ランダムにパウリ基底 (X, Y, Z) を選択して測定します。
    //
    for (int i = 0; i < number_of_total_measurements; i++) {
      for (int j = 0; j < system_size; j++) {
        printf("%c ", Pauli[rand() % 3]);
      }
      printf("\n");
    }
  }
  //
  // 古典シャドウの非ランダム化バージョンを実行
  //
  else if (strcmp(argv[1], "-d") == 0) {
    read_all_observables(argv[3]);

    //
    // 非ランダム化プロセスの効率的な使用のためにいくつかの定数を事前計算
    //
    double expm1eta = expm1(-eta / 2); // expm1eta = e^(-eta / 2) - 1
    for (int k = 0; k < max_k_local + 1; k++) {
      log1ppow1o3k.push_back(log1p(pow(1.0 / 3.0, k) * expm1eta));
    }

    //
    // 各局所観測量をこれだけの回数測定したい
    //
    number_of_measurements_per_observable = stoi(argv[2]);

    //
    // 古典シャドウの非ランダム化バージョン
    //
    // ランダムに選ぶ代わりに、未測定の観測量を効率的にカバーできるような
    // パウリ基底を決定論的に (貪欲法的に) 選択します。
    //

    // すべての観測量について、
    // 以前のすべての測定繰り返しの中で、その観測量が何回測定されたか
    vector<int> cur_num_of_measurements; // "現在の測定回数" を意味します
    cur_num_of_measurements.resize(number_of_observables, 0); // 0 で初期化

    // すべての観測量について、
    // 現在の測定繰り返しにおいて、その観測量を測定するために
    // いくつのパウリ演算子が一致する必要があるか
    vector<int> how_many_pauli_to_match;
    how_many_pauli_to_match.resize(number_of_observables);

    for (int measurement_repetition = 0; measurement_repetition < INF;
         measurement_repetition++) {
      for (int i = 0; i < (int)observables.size(); i++)
        how_many_pauli_to_match[i] =
            observables[i].size(); // k-local 観測量の場合は k で初期化

      double shift = (sum_cnt == 0) ? 0 : sum_log_value / sum_cnt;
      sum_log_value = 0.0;
      sum_cnt = 0;

      for (int ith_qubit = 0; ith_qubit < system_size; ith_qubit++) {
        double prob_of_failure[3]; // X, Y, または Z を選ぶための失敗確率
        double smallest_prob_of_failure = -1;

        //
        // 現在の繰り返しで ith_qubit に対してパウリ測定を選ぶ場合
        //
        for (int pauli = 0; pauli < 3; pauli++) {
          prob_of_failure[pauli] = 0;

          // すべてのパウリ観測量 p について、スコアを計算できる
          for (int p = 0; p < 3; p++) {
            for (int i : observables_acting_on_ith_qubit[ith_qubit][p]) {
              if (pauli == p) {
                int pauli_to_match_next_step =
                    how_many_pauli_to_match[i] == INF
                        ? INF
                        : how_many_pauli_to_match[i] - 1;
                double prob_next_step = fail_prob_pessimistic(
                    cur_num_of_measurements[i], pauli_to_match_next_step,
                    observables_weight[i], shift);
                double prob_current_step = fail_prob_pessimistic(
                    cur_num_of_measurements[i], how_many_pauli_to_match[i],
                    observables_weight[i], shift);
                prob_of_failure[pauli] += prob_next_step - prob_current_step;
              } else {
                double prob_next_step =
                    fail_prob_pessimistic(cur_num_of_measurements[i], INF,
                                          observables_weight[i], shift);
                double prob_current_step = fail_prob_pessimistic(
                    cur_num_of_measurements[i], how_many_pauli_to_match[i],
                    observables_weight[i], shift);
                prob_of_failure[pauli] += prob_next_step - prob_current_step;
              }
            }
          }

          if (smallest_prob_of_failure == -1)
            smallest_prob_of_failure = prob_of_failure[pauli];
          else
            smallest_prob_of_failure =
                min(smallest_prob_of_failure, prob_of_failure[pauli]);
        }

        // 最も低い失敗確率を持つものを選ぶ
        int the_best_pauli = 0;
        for (int pauli = 0; pauli < 3; pauli++) {
          if (smallest_prob_of_failure == prob_of_failure[pauli]) {
            printf("%c ", 'X' + pauli);
            the_best_pauli = pauli;
            break;
          }
        }

        for (int pauli = 0; pauli <= 2; pauli++) {
          for (int i : observables_acting_on_ith_qubit[ith_qubit][pauli]) {
            if (the_best_pauli == pauli) {
              if (how_many_pauli_to_match[i] != INF)
                how_many_pauli_to_match[i] -= 1;
            } else
              how_many_pauli_to_match[i] = INF;
          }
        }
      }
      printf("\n");

      for (int i = 0; i < (int)observables.size(); i++)
        if (how_many_pauli_to_match[i] == 0)
          cur_num_of_measurements[i]++;

      //
      // すべての観測量の測定回数をチェック
      //
      int success = 0;
      for (int i = 0; i < (int)observables.size(); i++)
        if (cur_num_of_measurements[i] >=
            floor(observables_weight[i] *
                  number_of_measurements_per_observable))
          success += 1;
      fprintf(stderr, "[Status %d: %d]\n", measurement_repetition + 1, success);

      if (success == (int)observables.size())
        break;
    }
  }
}
