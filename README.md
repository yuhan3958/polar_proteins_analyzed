## 전체 연구 원클릭 실행

모든 연구 분석을 한 번에 실행하려면 다음 명령어를 사용하세요.

```bash
python -m pip install -r requirements.txt
python run_all.py
```

**필수 외부 도구:** MAFFT, FastTree가 실행 경로(PATH)에 설치되어 있어야 합니다. 구조 분석에는 RCSB PDB에 접속할 인터넷이 필요합니다. 도구가 없으면 조용히 분석을 생략하지 않고 실패 원인을 출력합니다.

```bash
python run_all.py --skip-phylogeny --skip-structure
```

위 명령어는 외부 계통 분석 도구와 PDB 다운로드를 제외한 분석에 사용합니다. 순차 실행의 성공·실패와 결과 파일 경로는 `dataset_results/full_run_manifest.json`에 저장됩니다. 개별 스크립트를 새로 추가하면 `run_all.py`의 분석 단계와 CI도 반드시 갱신해야 합니다.

# GH2 단백질 비교 분석

### Windows: MAFFT/FastTree 오류

`python run_all.py`에서 `MAFFT and FastTree are required`가 나오면 **Python pip 패키지가 아닌 외부 실행 도구**가 PATH에서 보이지 않는 것입니다.

권장 방법은 WSL2 Ubuntu에서 두 도구와 Python 환경을 설치해 실행하는 것입니다.

```powershell
wsl --install -d Ubuntu
```

Ubuntu 터미널에서:

```bash
sudo apt update
sudo apt install -y mafft fasttree python3-venv
# Windows 프로젝트 폴더로 이동 (C: 경로를 /mnt/c/...로 변환)
cd /mnt/c/dev/PythonProject/PythonProject/polar_proteins_analyze
python3 -m venv .venv-wsl
source .venv-wsl/bin/activate
python -m pip install -r requirements.txt
python run_all.py
```

WSL 설치 후에는 Ubuntu를 다시 열어야 할 수 있습니다. Windows에서 내려받은 실행 파일을 사용한다면 `mafft`와 `FastTree` 명령이 Python 프로세스의 PATH에서 실제로 실행 가능한지도 확인하세요.

**도구가 없을 때도** 새 통합 실행기는 계통 분석만 미완료로 기록하고, 독립된 3D 구조 분석을 계속 수행합니다. `dataset_results/full_run_manifest.json`에 실패·차단·성공 상태가 남고 실행 종료 코드는 1이므로 부분 완료를 전체 성공으로 오인하지 않습니다.

## 실행

PyCharm에서 **`run_analysis.py`를 우클릭 → Run** 하세요. 추가 실행 인자는 필요 없습니다.
준비된 저온군 10개와 중온군 10개를 읽어 190쌍을 비교하고, 프로젝트의 `dataset_results/`에 결과를 저장합니다.
다시 실행하면 같은 이름의 결과 파일을 갱신합니다.

프로젝트 폴더의 PowerShell에서도 실행할 수 있습니다.

```powershell
.\.venv\Scripts\python.exe run_analysis.py
```

## 결과

| 파일 | 내용 |
|---|---|
| `dataset_results/features.csv` | 20개 단백질의 아미노산 조성, 분자량·pI·GRAVY 추정값 |
| `dataset_results/group_summary.csv` | 저온군·중온군 조성 비교와 부분집합 비교 |
| `dataset_results/pairwise_summary.csv` | 190쌍의 서열 일치율, 치환·갭 요약 |
| `dataset_results/pairwise/` | 각 쌍의 정렬 FASTA와 치환 위치 |
| `dataset_results/windows.csv` | 정렬의 50열 구간별 비교 |
| `dataset_results/composition_group.png` | 아미노산 20종의 군별 평균·표준편차 |
| `dataset_results/GPR_distribution.png` | G·P·R의 개별 단백질 분포 |
| `dataset_results/identity_heatmap.png` | 20개 단백질의 서열 일치율 행렬 |
| `dataset_results/summary.json` | 분석 설정, 입력 해시, 버전, 한계 |

## 스크립트 구분

- `run_analysis.py`: 준비된 20개 단백질 전체 분석.
- `analyze_polar_proteins.py`: 인자 없이 실행하면 기존 두 단백질 비교. 전체 분석은 `--metadata data/metadata.csv --out dataset_results` 지정.
- `scripts/research_gh2_sources.py`: 후보와 문헌 자료 수집.
- `scripts/collect_gh2_dataset.py`: 선정 목록에 따른 서열 다운로드·검증. `--offline`은 저장된 원본으로 입력을 재생성.

## 해석 범위

저온군에는 저온 내성균과 같은 속의 여러 균주가 포함됩니다. 일부 분류는 배양 조건이나 종 수준 근거에 의존하며, 배양 온도는 성장 최적 온도가 아닙니다. 결과는 선정된 표본의 기술적 비교이며, 계통 차이를 보정한 통계적 유의성이나 저온 적응의 인과관계를 입증하지 않습니다.

## 소논문 초안 자동 생성

`python run_all.py` 실행 시 모든 분석 단계가 성공하면 마지막으로 `paper/gh2_exploratory_manuscript_ko.md`가 자동 생성됩니다. 원고는 고정 숫자를 복사하지 않고 `dataset_results/`의 실제 CSV·JSON 결과를 읽습니다. 독립 검증 전 생물학적 인과관계나 유의성을 주장하지 않는 **탐색적 소논문**입니다. `python scripts/generate_manuscript.py`만 재실행할 수도 있습니다.

