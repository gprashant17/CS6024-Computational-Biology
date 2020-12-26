#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm> 
#include <functional>

using namespace std;
using std::vector;

string rtrim(const string&);
vector<string> split(const string&);
vector<int> GibbsSampler(string*, int, int, int);
vector<int> initializeMotifs(string*, int, int);
map<char, vector<float>> MotifToProfile(vector<int>, string*, int, int);
map<char, vector<float>> ProfileDelete(map<char, vector<float>>, string, int, int);
map<char, vector<float>> ProfileAdd(map<char, vector<float>>, string, int, int);
float CalculateProb(map<char, vector<float>>, string, int);
int InverseSampling(vector<float>, int, float);
int Score(map<char, vector<float>>, int, int);

int rand();
int rseed = 0;

#define RAND_MAX ((1U << 31) - 1)

inline int rand() {
    return rseed = (rseed * 1103515245 + 12345) & RAND_MAX;
}


vector<int> GibbsSampler(string* sequences, int k, int t, int N) {
    vector<int> Motifs(t, 0);
    map<char, vector<float>> Profile;

    vector<int> BestMotifs(Motifs);
    map<char, vector<float>> BestProfile;
    map<vector<int>, int> score_map;

    Motifs = initializeMotifs(sequences, k, t);
    Profile = MotifToProfile(Motifs, sequences, k, t);
    int score_profile = Score(Profile, k, t);

    int BestScore = score_profile;

    int d, n;

    for (int rstart = 0; rstart < 40; rstart++) {

        Motifs = initializeMotifs(sequences, k, t);
        Profile = MotifToProfile(Motifs, sequences, k, t);

        for (int iter = 0; iter < N; iter++) {

            d = rand() % t;
            map<char, vector<float>> ProfileDel = ProfileDelete(Profile, sequences[d].substr(Motifs[d], k), k, t);

            n = sequences[d].size() - k + 1;

            vector<float> Pr(n, 0);
            float sum = 0;
            for (int i = 0; i < n; i++) {
                Pr[i] = CalculateProb(ProfileDel, sequences[d].substr(i, k), k);
                sum += Pr[i];
            }

            int pos = InverseSampling(Pr, n, sum);

            Motifs[d] = pos;

            Profile = ProfileAdd(ProfileDel, sequences[d].substr(pos, k), k, t);

            if (score_map.find(Motifs) != score_map.end()) {
                score_profile = score_map[Motifs];
            }
            else {
                score_profile = Score(Profile, k, t);
            }

            if (score_profile < BestScore) {
                for (int m = 0; m < t; m++) {
                    BestMotifs[m] = Motifs[m];
                }
                BestScore = score_profile;
            }
        }
    }
    return BestMotifs;
}

vector<int> initializeMotifs(string* sequences, int k, int t) {
    vector<int> Motifs(t, 0);
    int n;
    for (int i = 0; i < t; i++) {
        n = sequences[i].size() - k + 1;
        Motifs[i] = rand() % n;
    }
    return Motifs;
}

map<char, vector<float>> MotifToProfile(vector<int> Motifs, string* sequences, int k, int t) {
    map<char, vector<float>> Profile;

    vector<float> vec(k, (1 / (float)(t + 4)));

    Profile['A'] = vec;
    Profile['T'] = vec;
    Profile['C'] = vec;
    Profile['G'] = vec;

    for (int i = 0; i < t; i++) {
        for (int j = 0; j < k; j++) {
            Profile[char(sequences[i].substr(Motifs[i], k)[j])][j] += 1 / (float)(t + 4);
        }
    }
    return Profile;
}

map<char, vector<float>> ProfileDelete(map<char, vector<float>> Profile, string pattern, int k, int t) {
    vector<char> bases = { 'A','T','C','G' };

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < 4; j++) {
            if (bases[j] != pattern[i]) {
                Profile[bases[j]][i] *= (t + 4) / (float)(t + 3);
            }
            else {
                Profile[bases[j]][i] = (Profile[bases[j]][i] * (t + 4) - 1) / (float)(t + 3);
            }
        }
    }
    return Profile;
}

map<char, vector<float>> ProfileAdd(map<char, vector<float>> Profile, string pattern, int k, int t) {
    vector<char> bases = { 'A','T','C','G' };

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < 4; j++) {
            if (bases[j] != pattern[i]) {
                Profile[bases[j]][i] *= (t + 3) / (float)(t + 4);
            }
            else {
                Profile[bases[j]][i] = (Profile[bases[j]][i] * (t + 3) + 1) / (float)(t + 4);
            }
        }
    }
    return Profile;
}

float CalculateProb(map<char, vector<float>> Profile, string pattern, int k) {
    float prob = 1;
    for (int i = 0; i < k; i++) {
        prob *= Profile[pattern[i]][i];
    }
    return prob;
}

int InverseSampling(vector<float> Pr, int len, float sum) {


    int* percentage;
    percentage = new int[len];
    int sum_per = 0;

    int distribution[100];
    int index = 0;
    for (int i = 0; i < len; i++)
    {
        percentage[i] = (Pr[i] / sum) * 100;
        sum_per += percentage[i];
        int times = percentage[i];
        while (times-- > 0)
        {
            distribution[index] = i;
            index++;
        }
    }
    //adjust the error at distribution[99]
    for (int i = sum_per; i < 100; i++)
    {

        distribution[i] = distribution[(rand() % sum_per)];
    }

    return distribution[rand() % 100];

}

int Score(map<char, vector<float>> Profile, int k, int t) {
    int Score = 0;

    vector<char> bases = { 'A','T','C','G' };

    for (int i = 0; i < k; i++) {
        int max = 0;
        for (int j = 0; j < 4; j++) {
            int s = Profile[bases[j]][i] * (4 + t) - 1;
            if (max < s) {
                max = s;
            }
        }
        Score += (t - max);
    }
    return Score;
}

int main() {
    // Write C++ code here
    string input_values;
    getline(cin, input_values);
    vector<string> input_values_final = split(rtrim(input_values));

    int k = stoi(input_values_final[0]);
    int t = stoi(input_values_final[1]);
    int N = stoi(input_values_final[2]);

    string* sequences;
    sequences = new string[t];

    for (int i = 0; i < t; i++) {
        getline(cin, sequences[i]);
    }
    vector<int> Motifs = GibbsSampler(sequences, k, t, N);

    for (int j = 0; j < t; j++) {
        cout << sequences[j].substr(Motifs[j], k) << "\n";
    }

    return 0;
}


string rtrim(const string& str) {
    string s(str);

    s.erase(
        find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(),
        s.end()
    );

    return s;
}

vector<string> split(const string& str) {
    vector<string> tokens;

    string::size_type start = 0;
    string::size_type end = 0;

    while ((end = str.find(" ", start)) != string::npos) {
        tokens.push_back(str.substr(start, end - start));

        start = end + 1;
    }

    tokens.push_back(str.substr(start));

    return tokens;
}   