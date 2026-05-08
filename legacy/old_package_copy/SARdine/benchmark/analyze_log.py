#!/usr/bin/env python3
"""Quick log analysis script."""
import re

# Read log
with open('/home/datacube/apps/SARdine/SARdine/benchmark/results/benchmark_log.txt', 'r') as f:
    lines = f.readlines()

# Parse step timings
step_times = {}
step_pattern = re.compile(r'STEP (\d+): (.+?) \(([\d.]+)s\)')

for line in lines:
    match = step_pattern.search(line)
    if match:
        step_num = int(match.group(1))
        step_name = match.group(2).strip()
        duration = float(match.group(3))
        step_times[step_num] = {
            'name': step_name,
            'duration_s': duration,
        }

# Print summary
print('Per-Stage Timings:')
print('-' * 60)
total = sum(st['duration_s'] for st in step_times.values())
for num in sorted(step_times.keys()):
    st = step_times[num]
    pct = st['duration_s'] / total * 100 if total > 0 else 0
    print(f"Step {num:2d}: {st['name']:40s} {st['duration_s']:7.1f}s ({pct:5.1f}%)")

print('-' * 60)
print(f'Total: {total:.1f}s')
print()
print('Top 5 Hotspots:')
sorted_steps = sorted(step_times.items(), key=lambda x: x[1]['duration_s'], reverse=True)[:5]
for i, (step_num, st) in enumerate(sorted_steps, 1):
    print(f"{i}. Step {step_num}: {st['name']} - {st['duration_s']:.1f}s ({st['duration_s']/total*100:.1f}%)")
