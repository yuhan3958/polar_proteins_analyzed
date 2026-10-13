"""Generate Korean research manuscript from validated analysis outputs."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "dataset_results"
PAPER = ROOT / "paper" / "gh2_exploratory_manuscript_ko.md"

def load_json(rel):
    return json.loads((DATA / rel).read_text(encoding="utf-8"))

def load_csv(rel):
    with (DATA / rel).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def make_paper():
    dataset = load_json("summary.json")
    contrasts = {r["amino_acid"]:r for r in load_csv("audit/group_sensitivity.csv")}
    phy = load_json("phylogeny/summary.json")
    supports = load_json("phylogeny/support_audit.json")
    structures = load_json("structures/summary.json")
    thermal = load_json("extended/summary.json")["group_metadata"]
    if len(contrasts) != 3 or dataset["total_pairs"] != 190:
        raise ValueError("Unexpected dataset version: manuscript requires 20 proteins")
    tab = []
    for aa in ("G","P","R"):
        r=contrasts[aa]
        tab.append(f"| {aa} | {float(r['protein_weighted_difference_pp']):+.3f} | {float(r['equal_genus_difference_pp']):+.3f} |")
    st = []
    for r in structures["structures"]:
        st.append(f"| {r['accession']} | {r['pdb']} | {r['mapped_residues']} | {100*r['fasta_to_structure_identity']:.2f}% |")
    text = f"""# 저온 적응 미생물과 중온성 미생물의 GH2 β-galactosidase 아미노산 조성 및 활성 부위 구조 비교

**형식:** 공개 생물정보학 자료에 기반한 탐색적 소논문 초안 (독립적 생물학적 실험 아님)

## 초록

저온 적응 효소의 아미노산 조성은 저온 활성과 열안정성에 관련된 후보 특성으로 논의되어 왔다. 본 연구에서는 공개 서열 중 저온군 {dataset['group_counts']['cold']}개와 중온군 {dataset['group_counts']['meso']}개의 GH2 β-galactosidase를 선정하고, {dataset['total_pairs']}쌍의 서열 정렬과 G(글리신)·P(프롤린)·R(아르기닌)의 조성 비교를 실시하였다. 단백질별 단순 평균과 미생물 속별 균형 평균을 비교하였으며, MAFFT 다중 서열 정렬, FastTree 단백질 계통 추정, 결정구조 2개의 촉매 관련 잔기 거리 계산을 수행하였다. 글리신의 저온군−중온군 평균 차이는 단순 평균 {float(contrasts['G']['protein_weighted_difference_pp']):+.3f}pp에서 속별 균형 후 {float(contrasts['G']['equal_genus_difference_pp']):+.3f}pp로 감소하였고 프롤린은 부호가 역전되었다. 계통 구성에 민감한 결과로 판단되며 저온 적응의 인과적 기작은 확정할 수 없다.

**주요어:** β-galactosidase, GH2, 저온 적응, 아미노산 조성, 계통 편향, 촉매 잔기

## 1. 서론

효소의 저온 적응은 낮은 온도에서의 반응 효율과 구조적 안정성 사이의 균형 문제를 제기한다. Hoyoux 등(2001)은 남극 *Pseudoalteromonas haloplanktis* TAE79 β-galactosidase가 대장균 효소보다 낮은 온도에서 높은 효율을 보이는 사례를 보고하고, 아르기닌·프롤린 조성과 촉매 부위 주변 글리신의 역할을 논의하였다[1]. 그러나 서로 먼 계통에 속한 소수의 효소에서 관찰한 패턴이 저온성 미생물 전체에 재현되는지는 별도 확인이 필요하다.

연구 질문은 (1) 저온군에서 G가 많고 P와 R이 적은가, (2) 이 차이는 속별 표본 불균형에 민감한가, (3) 촉매 관련 잔기의 상동 위치가 실험 구조에서 측정 가능한가이다.

## 2. 재료 및 방법

### 2.1 단백질 및 메타데이터

데이터셋은 저온군 {dataset['group_counts']['cold']}개, 중온군 {dataset['group_counts']['meso']}개 단백질로 구성했다. 접근번호, 계열 정보, 서열 해시 및 분류 출처는 저장소의 data/metadata.csv를 사용했다. 효소 활성 최적 온도, 배양 온도, 미생물 성장 최적 온도는 서로 다른 물리량으로 취급하였다. 성장 최적 온도 수치가 기록된 표본은 저온군 {thermal['cold']['growth_optimum_recorded']}개, 중온군 {thermal['meso']['growth_optimum_recorded']}개였다. URL의 존재는 원문 수치 검증 완료를 의미하지 않는다.

### 2.2 서열 조성 및 계통 민감도

20개 표준 아미노산의 빈도를 각 단백질 길이로 나눈 비율을 계산하였다. G/P/R의 집단별 평균과 함께, 속마다 하나의 평균을 먼저 구한 뒤 속에 동일한 가중치를 주는 민감도 분석을 수행하였다. 이는 계통학적 비교모형이 아니라 과대표집에 관한 기술적 진단이다.

### 2.3 다중 정렬·계통 추정

MAFFT {phy['tools']['mafft'][0]}를 이용해 총 {phy['alignment_columns']}열의 정렬을 얻고 FastTree -wag으로 근사 최대우도 단백질 계통수를 추정하였다. FastTree SH-like 국소 지지도는 일반적인 부트스트랩 반복 지지도와 구별하였다.

### 2.4 활성 관련 위치 및 입체구조

대장균 LacZ의 Arg388, Glu461, Met502, Tyr503, Glu537을 기준으로 상동 잔기를 서열 정렬로 찾았다. 실험 결정구조 1DP0(대장균)과 1YQ2(*Arthrobacter* C2-2)에서 Cα–Cα 거리 및 가장 가까운 무거운 원자 거리를 산출하였다[2]. 정적 구조의 거리 차이는 효소 유연성의 직접 측정이 아니다.

## 3. 결과

### 3.1 아미노산 조성

표 1. 저온군 평균에서 중온군 평균을 뺀 차이(퍼센트포인트).

| 잔기 | 단백질별 평균 차이 | 속별 동일 가중 차이 |
|---|---:|---:|
{chr(10).join(tab)}

저온군에서는 *Arthrobacter*가 6/10개를 차지하는 반면 중온군은 10개가 서로 다른 속에 속한다. 따라서 G 차이가 속 균형 조정 후 거의 사라지고 P의 부호가 바뀐 점은 표본 구성에 대한 민감성을 보여 준다.

### 3.2 계통 추정 및 지지도

단백질 {supports['leaf_count']}개의 추정 나무에는 내부 노드 {supports['internal_nodes']}개 중 국소 지지도 기록이 있는 노드 {supports['nodes_with_support']}개가 포함된다. 지지도 0.70 미만은 {supports['local_support_below_0_70']}개, 0.90 미만은 {supports['local_support_below_0_90']}개이며 최소값은 {supports['local_support_min']:.2f}이다. 이 나무는 종 계통수가 아닌 하나의 효소 계통수다.

### 3.3 결정구조

표 2. 결정구조에 매핑된 촉매 관련 위치.

| 단백질 | PDB | 대응 잔기 수 | FASTA–구조 서열 일치율 |
|---|---|---:|---:|
{chr(10).join(st)}

총 {structures['n_distances']}개의 잔기 쌍 거리 측정값은 dataset_results/structures/distances.csv에 수록하였다. 물리적 거리는 서로 다른 효소의 구조적 환경 비교에 활용할 수 있으나, 여기서는 계통·올리고머화·결정 조건 차이가 함께 존재한다.

## 4. 논의

원래 가설은 저온 적응 효소의 G 증가와 P/R 감소였으나, 실제 결과는 표본 균형 방식에 크게 의존하였다. 특히 글리신의 관찰된 증가분은 속별 균형 조정에서 거의 없어졌다. 이러한 민감성은 특정 아미노산 조성을 온도 적응의 보편적 지표로 취급하는 데 주의를 요구한다. *Arthrobacter* C2-2의 β-galactosidase는 대장균 효소와 올리고머 구조가 달라 단순한 국소 거리 비교에도 제한이 있다[2].

## 5. 한계

본 연구는 비무작위 표본 20개에 대한 탐색적 분석이다. 미생물 성장 최적 온도가 상당수 누락되어 온도를 연속 변수로 회귀할 수 없고, 같은 속의 표본이 독립적이라는 보장이 없다. 서열별 평균 차이에서 p값이나 인과 효과를 주장하지 않으며, 계통수의 SH-like 지지도는 일반적인 부트스트랩 신뢰구간이 아니다. 활성 부위 관련 잔기의 3차원 거리는 측정하였으나 실제 저온 반응속도나 구조 동역학은 측정하지 않았다.

## 6. 결론

현재 검증된 주요 결론은 **선정된 GH2 β-galactosidase 표본의 G/P/R 평균 차이가 속별 표본 구성에 민감하다**는 것이다. 따라서 저온 적응의 보편적 분자 기작은 검증 전 가설로 남겨야 한다.

## 참고문헌

[1] Hoyoux, A. et al. (2001). Cold-Adapted β-Galactosidase from the Antarctic Psychrophile *Pseudoalteromonas haloplanktis*. *Applied and Environmental Microbiology*, 67, 1529–1535. DOI: 10.1128/AEM.67.4.1529-1535.2001.

[2] Skálová, T. et al. (2005). Cold-active beta-galactosidase from *Arthrobacter* sp. C2-2 forms compact 660 kDa hexamers. *Journal of Molecular Biology*, 353, 282–294. DOI: 10.1016/j.jmb.2005.08.028.

## 재현성 및 데이터 공개

연구 코드와 표본 선정 메타데이터는 https://github.com/yuhan3958/polar_proteins_analyze 에 정리되어 있다. 모든 파이프라인은 python run_all.py로 실행하며, 본 원고는 python scripts/generate_manuscript.py로 동일한 분석 결과 파일에서 재생성한다.

**원고 상태:** 소논문 형태의 초안. 개별 표본의 모든 온도 출처와 해석은 최종 제출 전 독립 검토가 필요하다.
"""
    PAPER.parent.mkdir(parents=True,exist_ok=True)
    PAPER.write_text(text,encoding="utf-8")
    print("Wrote", PAPER)

if __name__=="__main__":
    make_paper()
