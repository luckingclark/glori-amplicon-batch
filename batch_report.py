"""Validate GLORI counts, left-join reference A positions, and render offline reports."""
import csv
import html
import math
from pathlib import Path
from urllib.parse import quote

from batch_common import read_json, write_tsv, now

FDR_FIELDS = ['Chr', 'Sites', 'Strand', 'Gene', 'CR', 'AGcov', 'Acov', 'Genecov',
              'Ratio', 'Pvalue', 'P_adjust']
SITE_FIELDS = ['Sample', 'Chr', 'Sites', 'Strand', 'Gene', 'Cutoff', 'Total_all',
               'AGcov_all', 'Acov_all', 'Gcov_all', 'Ratio_all', 'Percent_all',
               'Total_used', 'AGcov', 'Acov', 'Gcov', 'Ratio', 'Percent', 'Signal',
               'Background_CR', 'FDR_status', 'Pvalue', 'P_adjust', 'Observation_status']
QC_FIELDS = ['Sample', 'Status', 'Raw_reads', 'Clean_reads', 'Retention_pct',
             'Raw_Q30_pct', 'Clean_Q30_pct', 'Mapped_reads', 'Mapping_pct',
             'Reverse_reads', 'A_sites', 'Observed_A_sites', 'Significant_sites',
             'Raw_FastQC', 'Clean_FastQC', 'Seconds', 'Warnings', 'Error', 'Attempt']


def ratio(a, ag, scale=1):
    return scale * a / ag if ag else 'NA'


def counts(values):
    total, ag, a = map(int, values)
    if not 0 <= a <= ag <= total:
        raise ValueError(f'无效计数: total={total}, AG={ag}, A={a}')
    return total, ag, a


def probability(value):
    number = float(value)
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise ValueError(f'无效比例/P 值: {value}')
    return number


def site_key(row):
    return row['Chr'], int(row['Sites']), row['Strand']


def reference_keys(records):
    return {(name, position, '+') for name, sequence in records.items()
            for position, base in enumerate(sequence, 1) if base == 'A'}


def read_fdr(path, allowed):
    result = {}
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != FDR_FIELDS:
            raise ValueError(f'FDR 表头不完整: {path}')
        for row in reader:
            key = site_key(row)
            if key not in allowed or key in result or None in row or None in row.values():
                raise ValueError(f'FDR 位点重复或与参考不匹配: {key}')
            ag, a = int(row['AGcov']), int(row['Acov'])
            if not 0 <= a <= ag or ag == 0:
                raise ValueError(f'FDR 计数不合法: {key}')
            for field in ['CR', 'Ratio', 'Pvalue', 'P_adjust']:
                probability(row[field])
            if abs(float(row['Ratio']) - a / ag) > 0.000011:
                raise ValueError(f'FDR 比例与计数不一致: {key}')
            result[key] = row
    return result


def export_sites(sample, output, destination, cutoff='3'):
    """No re-testing: significant calls retain the upstream per-sample FDR."""
    output, destination = Path(output), Path(destination)
    prefix = sample['sample']
    allowed = reference_keys(sample['records'])
    significant = read_fdr(output / f'{prefix}.totalm6A.FDR.csv', allowed)
    backgrounds = {}
    with (output / f'{prefix}.totalCR.txt').open(encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['SA', 'A-to-G_ratio']:
            raise ValueError('转换率表头不正确。')
        for row in reader:
            if row['SA'] in {'#' + key for key in sample['records']}:
                name = row['SA'][1:]
                if name in backgrounds:
                    raise ValueError(f'重复的转换率: {name}')
                backgrounds[name] = probability(row['A-to-G_ratio'])
    observed = {}
    path = output / 'tmp' / f'{prefix}.totalformat.txt'
    with path.open(encoding='utf-8') as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if not row or row[0].startswith('#'):
                continue
            if len(row) != 26:
                raise ValueError(f'位点计数文件列数错误: {path}')
            key = row[0], int(row[1]), row[2]
            if key not in allowed or key in observed:
                raise ValueError(f'计数位点重复或与参考不匹配: {key}')
            total, ag, a = counts(row[8:11])
            if int(row[13]) != ag - a:
                raise ValueError(f'G 计数不一致: {key}')
            bins = dict(field.split(';', 1) for field in row[14:26])
            used_total, used_ag, used_a = counts(bins[cutoff].split(','))
            if (used_total > total or used_ag > ag or used_a > a
                    or used_ag - used_a > ag - a):
                raise ValueError(f'过滤后的计数大于过滤前: {key}')
            observed[key] = dict(Total_all=total, AGcov_all=ag, Acov_all=a, Gcov_all=ag-a,
                                 Ratio_all=ratio(a, ag), Percent_all=ratio(a, ag, 100),
                                 Total_used=used_total, AGcov=used_ag, Acov=used_a,
                                 Gcov=used_ag-used_a, Ratio=ratio(used_a, used_ag),
                                 Percent=ratio(used_a, used_ag, 100),
                                 Signal=ratio(used_total, total),
                                 Observation_status='observed' if used_ag else 'no_AG_after_cutoff')
    for key, call in significant.items():
        data = observed.get(key)
        if data is None or (data['AGcov'], data['Acov']) != (int(call['AGcov']), int(call['Acov'])):
            raise ValueError(f'FDR 位点与原始计数表不一致: {key}')
    rows = []
    for row in empty_sites(sample, 'no_usable_record'):
        key = site_key(row)
        row['Background_CR'] = backgrounds.get(row['Chr'], 'NA')
        if key in observed:
            row.update(observed[key])
            row['FDR_status'] = 'not_reported'
        if key in significant:
            row.update(FDR_status='significant', Pvalue=significant[key]['Pvalue'],
                       P_adjust=significant[key]['P_adjust'])
        rows.append(row)
    write_tsv(destination / 'all_A_sites.tsv', SITE_FIELDS, rows)
    write_tsv(destination / 'significant_sites.tsv', ['Sample'] + FDR_FIELDS,
              [dict(Sample=prefix, **row) for row in significant.values()])
    return dict(a_sites=len(allowed), observed_a_sites=len(observed), significant_sites=len(significant))


def empty_sites(sample, reason):
    for name, sequence in sample['records'].items():
        for position, base in enumerate(sequence, 1):
            if base != 'A':
                continue
            row = dict.fromkeys(SITE_FIELDS, 'NA')
            row.update(Sample=sample['sample'], Chr=name, Sites=position, Strand='+',
                       Gene=name, Cutoff='3', Observation_status=reason)
            yield row


def read_table(path, fields):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != fields:
            raise ValueError(f'结果表头不正确: {path}')
        rows = list(reader)
        if any(None in row or None in row.values() for row in rows):
            raise ValueError(f'结果表损坏: {path}')
        return rows


def collect(root, batch, scheduler_states=None, final=True):
    root = Path(root)
    scheduler_states = scheduler_states or {}
    qc_rows, all_rows, significant_rows = [], [], []
    links = []
    failures = 0
    for sample in batch['samples']:
        name = sample['sample']
        relative = batch['latest'][name]
        attempt = root / relative
        status_file = attempt / 'status.json'
        try:
            status = read_json(status_file) if status_file.exists() else {}
        except (ValueError, OSError) as error:
            status = dict(state='FAILED', error=f'状态文件损坏: {error}')
        state = status.get('state', 'NOT_STARTED')
        error = status.get('error', '')
        qc = status.get('qc', {})
        if state == 'SUCCESS':
            try:
                rows = read_table(attempt / 'all_A_sites.tsv', SITE_FIELDS)
                sig = read_table(attempt / 'significant_sites.tsv', ['Sample'] + FDR_FIELDS)
                if (len(rows) != len(reference_keys(sample['records']))
                        or {site_key(row) for row in rows} != reference_keys(sample['records'])
                        or any(row['Sample'] != name for row in rows + sig)):
                    raise ValueError('成功样品的汇总位点与参考或样品名不匹配。')
                all_rows.extend(rows)
                significant_rows.extend(sig)
            except (OSError, ValueError, KeyError) as problem:
                state, error = 'FAILED', f'结果文件缺失或损坏: {problem}'
        if state != 'SUCCESS':
            failures += 1
            if not error and final:
                state = 'INCOMPLETE'
                error = ('作业未写入成功状态；Slurm: ' + scheduler_states.get(name, '未知')
                         + '。请检查作业日志（可能超时、取消、内存不足或尚未运行）。')
            all_rows.extend(empty_sites(sample, 'sample_failed' if final or state == 'FAILED' else 'sample_pending'))
        row = dict.fromkeys(QC_FIELDS, 'NA')
        raw, clean = qc.get('raw', {}), qc.get('clean', {})
        row.update(Sample=name, Status=state, Raw_reads=raw.get('reads', 'NA'),
                   Clean_reads=clean.get('reads', 'NA'), Retention_pct=qc.get('retention_pct', 'NA'),
                   Raw_Q30_pct=raw.get('q30_pct', 'NA'), Clean_Q30_pct=clean.get('q30_pct', 'NA'),
                   Mapped_reads=qc.get('mapped_reads', 'NA'), Mapping_pct=qc.get('mapping_pct', 'NA'),
                   Reverse_reads=qc.get('reverse_reads', 'NA'), A_sites=len(reference_keys(sample['records'])),
                   Observed_A_sites=qc.get('observed_a_sites', 'NA') if state == 'SUCCESS' else 'NA',
                   Significant_sites=qc.get('significant_sites', 'NA') if state == 'SUCCESS' else 'NA',
                   Raw_FastQC=qc.get('raw_fastqc', 'NA'), Clean_FastQC=qc.get('clean_fastqc', 'NA'),
                   Seconds=status.get('seconds', 'NA'), Warnings='; '.join(qc.get('warnings', [])),
                   Error=error, Attempt=relative)
        qc_rows.append(row)
        for file in sorted(attempt.glob('qc/*/*_fastqc.html')):
            links.append((name + ' / ' + file.parent.name + ' FastQC', file.relative_to(root).as_posix()))
        for file in ['logs/pipeline.log', 'logs/commands.json', 'status.json', 'all_A_sites.tsv']:
            if (attempt / file).exists():
                links.append((name + ' / ' + file, (Path(relative) / file).as_posix()))
    write_tsv(root / 'qc_summary.tsv', QC_FIELDS, qc_rows)
    write_tsv(root / 'all_A_sites.tsv', SITE_FIELDS, all_rows)
    write_tsv(root / 'significant_sites.tsv', ['Sample'] + FDR_FIELDS, significant_rows)
    render_report(root, qc_rows, significant_rows, links, batch)
    return failures


def render_report(root, qc, significant, links, batch):
    escape = lambda text: html.escape(str(text), quote=True)
    def table(rows, fields):
        head = ''.join('<th>' + escape(field) + '</th>' for field in fields)
        body = ''.join('<tr>' + ''.join('<td>' + escape(row.get(field, 'NA')) + '</td>'
                                      for field in fields) + '</tr>' for row in rows)
        return '<div class="scroll"><table><thead><tr>' + head + '</tr></thead><tbody>' + body + '</tbody></table></div>'
    downloads = [(name, name) for name in ['qc_summary.tsv', 'all_A_sites.tsv', 'significant_sites.tsv']]
    anchors = ''.join('<li><a href="' + quote(path, safe='/') + '">' + escape(label) + '</a></li>'
                      for label, path in downloads + links)
    success = sum(row['Status'] == 'SUCCESS' for row in qc)
    content = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>GLORI 批次报告</title>
<style>body{font:16px/1.65 system-ui,sans-serif;margin:32px auto;max-width:1280px;padding:0 20px;color:#203041;background:#f7f9fc}h1,h2{color:#123b58}.card{background:white;border:1px solid #dce4ec;padding:20px;border-radius:10px;margin:18px 0}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border-bottom:1px solid #ddd;padding:9px;text-align:left;white-space:nowrap}th{background:#eaf2f8}a{color:#075e96}code{background:#eef2f6;padding:2px 5px}</style>
<h1>GLORI 扩增子 m6A 分析报告</h1>'''
    content += f'<p>生成时间（UTC）：{escape(now())}。成功 {success} / {len(qc)} 个样品。</p>'
    content += '<div class="card"><b>如何读结果</b><p>Ratio = A / (A + G)，Percent 为百分比。'
    content += 'NA 表示无可用统计值；不是 0% 修饰。比例受比对、碱基质量、A-cutoff 和最大深度限制，未做背景或 PCR 校正。</p>'
    content += '<p>全部 A 位点请查看 all_A_sites.tsv。significant 表示通过原 GLORI 筛选；not_reported 表示未进入最终显著表，可能未满足前置筛选，不能等同于“没有修饰”。各样品独立进行 FDR。</p>'
    content += '<p>扩增子常有较高重复率和序列偏好；FastQC 的 WARN/FAIL 是 QC 提示，不自动判定实验失败，也不触发去重复。</p></div>'
    content += '<h2>样品 QC 与运行状态</h2>' + table(qc, ['Sample', 'Status', 'Raw_reads', 'Clean_reads',
                   'Retention_pct', 'Mapping_pct', 'Reverse_reads', 'Significant_sites', 'Warnings', 'Error'])
    content += '<h2>显著位点（最多预览 200 行）</h2>'
    content += table(significant[:200], ['Sample', 'Chr', 'Sites', 'AGcov', 'Acov', 'Ratio', 'P_adjust'])
    content += '<h2>结果及详细 QC 文件</h2><ul>' + anchors + '</ul>'
    content += '<p>离线查看时请下载整个批次目录，以保留 FastQC 和日志的相对链接。资源设置：' + escape(batch['resources']) + '</p></html>'
    path = root / 'report.html'
    temporary = root / 'report.html.tmp'
    temporary.write_text(content, encoding='utf-8')
    temporary.replace(path)
