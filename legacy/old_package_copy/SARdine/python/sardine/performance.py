"""
SARdine Performance Monitoring and Optimization

Utilities for monitoring and optimizing parallel processing performance.
"""

import time
import os
import psutil
import threading
from dataclasses import dataclass
from typing import Dict, Any, List, Optional


@dataclass
class PerformanceMetrics:
    """Performance metrics for processing operations"""
    start_time: float
    end_time: float
    duration: float
    cpu_count: int
    threads_used: int
    memory_usage_mb: float
    peak_memory_mb: float
    chunk_size: int
    parallel_enabled: bool
    step_name: str
    data_size_mb: float
    throughput_mbps: float
    # Detailed timing breakdowns (optional)
    timing_breakdown: Optional[Dict[str, float]] = None


class PerformanceMonitor:
    """Monitor performance of parallel processing operations with detailed timing breakdowns"""
    
    def __init__(self, enable_monitoring: bool = True):
        self.enable_monitoring = enable_monitoring
        self.metrics: List[PerformanceMetrics] = []
        self.current_step = None
        self.start_time = None
        self.peak_memory = 0.0
        self._monitoring_thread = None
        self._stop_monitoring = False
        # Detailed timing breakdowns for current step
        self._timing_breakdown: Dict[str, float] = {}
        self._timing_checkpoints: Dict[str, float] = {}
        
    def start_step(self, step_name: str, data_size_mb: float = 0.0, 
                   parallel_enabled: bool = True, threads_used: int = 0, 
                   chunk_size: int = 1024):
        """Start monitoring a processing step"""
        if not self.enable_monitoring:
            return
            
        self.current_step = {
            'name': step_name,
            'start_time': time.time(),
            'data_size_mb': data_size_mb,
            'parallel_enabled': parallel_enabled,
            'threads_used': threads_used or os.cpu_count(),
            'chunk_size': chunk_size,
            'start_memory': self._get_memory_usage()
        }
        
        self.peak_memory = self.current_step['start_memory']
        
        # Start memory monitoring thread
        if self._monitoring_thread is None or not self._monitoring_thread.is_alive():
            self._stop_monitoring = False
            self._monitoring_thread = threading.Thread(target=self._monitor_memory)
            self._monitoring_thread.daemon = True
            self._monitoring_thread.start()
        
        print(f"🔍 Monitoring: {step_name}")
        if parallel_enabled:
            print(f"   🚀 Parallel mode: {threads_used} threads, chunk size: {chunk_size}")
        else:
            print(f"   🔄 Sequential mode")
        
    def end_step(self):
        """End monitoring current step"""
        if not self.enable_monitoring or not self.current_step:
            return
            
        end_time = time.time()
        duration = end_time - self.current_step['start_time']
        current_memory = self._get_memory_usage()
        
        # Calculate throughput
        throughput = (self.current_step['data_size_mb'] / duration) if duration > 0 else 0
        
        metrics = PerformanceMetrics(
            start_time=self.current_step['start_time'],
            end_time=end_time,
            duration=duration,
            cpu_count=os.cpu_count(),
            threads_used=self.current_step['threads_used'],
            memory_usage_mb=current_memory,
            peak_memory_mb=self.peak_memory,
            chunk_size=self.current_step['chunk_size'],
            parallel_enabled=self.current_step['parallel_enabled'],
            step_name=self.current_step['name'],
            data_size_mb=self.current_step['data_size_mb'],
            throughput_mbps=throughput,
            timing_breakdown=self._timing_breakdown.copy() if self._timing_breakdown else None
        )
        
        self.metrics.append(metrics)
        
        print(f"   ✅ Completed in {duration:.2f}s")
        if self.current_step['data_size_mb'] > 0:
            print(f"   📊 Throughput: {throughput:.1f} MB/s")
        print(f"   💾 Memory: {current_memory:.1f} MB (peak: {self.peak_memory:.1f} MB)")
        
        # Print timing breakdown if available
        if self._timing_breakdown:
            print(f"   ⏱️  Timing breakdown:")
            for phase, phase_duration in sorted(self._timing_breakdown.items(), key=lambda x: x[1], reverse=True):
                percentage = (phase_duration / duration * 100) if duration > 0 else 0
                print(f"      • {phase}: {phase_duration:.3f}s ({percentage:.1f}%)")
        
        # Reset for next step
        self.current_step = None
        self._timing_breakdown = {}
        self._timing_checkpoints = {}
        
    def stop_monitoring(self):
        """Stop all monitoring"""
        self._stop_monitoring = True
        if self._monitoring_thread and self._monitoring_thread.is_alive():
            self._monitoring_thread.join(timeout=1.0)
    
    def _monitor_memory(self):
        """Monitor memory usage in background thread"""
        while not self._stop_monitoring:
            try:
                current_memory = self._get_memory_usage()
                if current_memory > self.peak_memory:
                    self.peak_memory = current_memory
                time.sleep(0.1)  # Check every 100ms
            except Exception:
                break
    
    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB"""
        try:
            process = psutil.Process()
            return process.memory_info().rss / 1024 / 1024
        except Exception:
            return 0.0
    
    def record_timing_checkpoint(self, checkpoint_name: str):
        """Record a timing checkpoint for detailed breakdown analysis
        
        Usage:
            monitor.record_timing_checkpoint("lut_precomputation")
            # ... do work ...
            monitor.record_timing_checkpoint("hot_loop")
            # ... do work ...
            monitor.end_step()  # Will show breakdown: lut_precomputation, hot_loop, etc.
        """
        if not self.enable_monitoring or not self.current_step:
            return
        
        current_time = time.time()
        self._timing_checkpoints[checkpoint_name] = current_time
        
        # Calculate duration since last checkpoint
        if len(self._timing_checkpoints) > 1:
            # Get previous checkpoint
            sorted_checkpoints = sorted(self._timing_checkpoints.items(), key=lambda x: x[1])
            if len(sorted_checkpoints) >= 2:
                prev_name, prev_time = sorted_checkpoints[-2]
                phase_duration = current_time - prev_time
                phase_name = f"{prev_name} → {checkpoint_name}"
                self._timing_breakdown[phase_name] = phase_duration
        else:
            # First checkpoint - time since step start
            step_start = self.current_step['start_time']
            phase_duration = current_time - step_start
            self._timing_breakdown[f"start → {checkpoint_name}"] = phase_duration
    
    def record_phase_duration(self, phase_name: str, duration: float):
        """Directly record a phase duration (for phases measured externally)
        
        Args:
            phase_name: Name of the phase (e.g., "SLC_reading", "calibration_LUT")
            duration: Duration in seconds
        """
        if not self.enable_monitoring:
            return
        self._timing_breakdown[phase_name] = duration
    
    def get_summary(self) -> Dict[str, Any]:
        """Get performance summary"""
        if not self.metrics:
            return {"status": "no_metrics"}
        
        total_duration = sum(m.duration for m in self.metrics)
        parallel_steps = [m for m in self.metrics if m.parallel_enabled]
        sequential_steps = [m for m in self.metrics if not m.parallel_enabled]
        
        summary = {
            "total_steps": len(self.metrics),
            "total_duration_seconds": total_duration,
            "parallel_steps": len(parallel_steps),
            "sequential_steps": len(sequential_steps),
            "average_throughput_mbps": sum(m.throughput_mbps for m in self.metrics if m.throughput_mbps > 0) / len([m for m in self.metrics if m.throughput_mbps > 0]) if any(m.throughput_mbps > 0 for m in self.metrics) else 0,
            "peak_memory_mb": max(m.peak_memory_mb for m in self.metrics),
            "cpu_cores_available": os.cpu_count(),
            "steps": []
        }
        
        for m in self.metrics:
            summary["steps"].append({
                "name": m.step_name,
                "duration_seconds": m.duration,
                "parallel_enabled": m.parallel_enabled,
                "threads_used": m.threads_used,
                "throughput_mbps": m.throughput_mbps,
                "memory_mb": m.memory_usage_mb
            })
        
        return summary
    
    def print_summary(self):
        """Print performance summary"""
        summary = self.get_summary()
        
        if summary.get("status") == "no_metrics":
            print("📊 No performance metrics collected")
            return
        
        print("\n" + "="*80)
        print("📊 PERFORMANCE SUMMARY")
        print("="*80)
        
        print(f"🕐 Total processing time: {summary['total_duration_seconds']:.1f} seconds")
        print(f"🏃 Parallel steps: {summary['parallel_steps']}/{summary['total_steps']}")
        print(f"💻 CPU cores available: {summary['cpu_cores_available']}")
        print(f"💾 Peak memory usage: {summary['peak_memory_mb']:.1f} MB")
        
        if summary['average_throughput_mbps'] > 0:
            print(f"📊 Average throughput: {summary['average_throughput_mbps']:.1f} MB/s")
        
        print(f"\n📋 Step-by-step performance:")
        for step in summary['steps']:
            mode = "🚀 PARALLEL" if step['parallel_enabled'] else "🔄 SEQUENTIAL"
            threads = f" ({step['threads_used']} threads)" if step['parallel_enabled'] else ""
            throughput = f", {step['throughput_mbps']:.1f} MB/s" if step['throughput_mbps'] > 0 else ""
            print(f"   {step['name']}: {step['duration_seconds']:.1f}s {mode}{threads}{throughput}")
        
        # Performance recommendations
        self._print_recommendations(summary)
        
    def _print_recommendations(self, summary: Dict[str, Any]):
        """Print performance optimization recommendations"""
        print(f"\n💡 PERFORMANCE RECOMMENDATIONS:")
        
        parallel_ratio = summary['parallel_steps'] / summary['total_steps']
        if parallel_ratio < 0.8:
            print("   • Consider enabling parallel processing for more steps")
        
        if summary['peak_memory_mb'] > 8000:  # > 8GB
            print("   • High memory usage detected - consider reducing chunk size")
        elif summary['peak_memory_mb'] < 2000:  # < 2GB  
            print("   • Low memory usage - consider increasing chunk size for better performance")
        
        if summary['average_throughput_mbps'] < 50:
            print("   • Low throughput detected - check if all CPU cores are being utilized")
        
        print("   • For best performance: ensure sufficient RAM and use SSD storage")
        print("   • Use --num-threads to explicitly set thread count if auto-detection is suboptimal")


def optimize_chunk_size(data_size_mb: float, num_threads: int, 
                       memory_limit_mb: float = 8000) -> int:
    """Optimize chunk size based on data size, threads, and memory constraints"""
    
    # Base chunk size calculation
    base_chunk = max(64, int(data_size_mb / num_threads / 4))
    
    # Memory constraint
    max_chunk_for_memory = int(memory_limit_mb / num_threads / 8)  # Conservative estimate
    
    # Choose the smaller of the two
    optimal_chunk = min(base_chunk, max_chunk_for_memory)
    
    # Ensure reasonable bounds
    optimal_chunk = max(64, min(optimal_chunk, 4096))
    
    return optimal_chunk


def get_system_info() -> Dict[str, Any]:
    """Get system information for performance tuning"""
    try:
        cpu_count = os.cpu_count()
        memory = psutil.virtual_memory()
        
        return {
            "cpu_cores": cpu_count,
            "memory_total_gb": memory.total / 1024**3,
            "memory_available_gb": memory.available / 1024**3,
            "memory_percent_used": memory.percent,
            "recommended_threads": cpu_count,
            "recommended_chunk_size": optimize_chunk_size(1000, cpu_count),  # For 1GB typical data
        }
    except Exception as e:
        return {"error": str(e)}


def print_system_info():
    """Print system information for performance optimization"""
    info = get_system_info()
    
    print("💻 SYSTEM PERFORMANCE INFO")
    print("="*50)
    print(f"🏭 CPU cores: {info.get('cpu_cores', 'unknown')}")
    print(f"💾 Total RAM: {info.get('memory_total_gb', 0):.1f} GB")
    print(f"💾 Available RAM: {info.get('memory_available_gb', 0):.1f} GB")
    print(f"📊 Memory usage: {info.get('memory_percent_used', 0):.1f}%")
    print(f"🚀 Recommended threads: {info.get('recommended_threads', 'auto')}")
    print(f"📦 Recommended chunk size: {info.get('recommended_chunk_size', 1024)}")
    
    # Performance tips
    memory_gb = info.get('memory_total_gb', 0)
    if memory_gb < 8:
        print("\n⚠️  LOW MEMORY: Consider using smaller chunk sizes and sequential processing for large files")
    elif memory_gb > 32:
        print("\n🚀 HIGH MEMORY: You can use large chunk sizes and maximum parallelization")
    else:
        print("\n✅ ADEQUATE MEMORY: Default settings should work well")
