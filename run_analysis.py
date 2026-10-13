"""PyCharm에서 이 파일을 실행하면 준비된 20개 GH2 단백질을 모두 비교한다.

결과: 프로젝트 루트의 dataset_results 폴더.
데이터 수집용 scripts/research_gh2_sources.py와는 별도의 분석 진입점이다.
"""

from pathlib import Path

from dataset_analysis import run_dataset_analysis


def main():
    root = Path(__file__).resolve().parent
    metadata = root / 'data' / 'metadata.csv'
    if not metadata.exists():
        raise FileNotFoundError('data/metadata.csv가 없습니다. 먼저 scripts/collect_gh2_dataset.py로 데이터를 수집하세요.')
    run_dataset_analysis(metadata, root / 'dataset_results', window_size=50)
    print(f'\n분석 완료. 결과 폴더: {root / "dataset_results"}')


if __name__ == '__main__':
    main()
