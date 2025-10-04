# Bug Report Template

## Bug Summary
**Bug Title**: [Concise description of the bug]
**Severity**: [Critical/High/Medium/Low]
**Priority**: [Urgent/High/Medium/Low]
**Status**: [New/Investigating/In Progress/Fixed/Closed]
**Assigned To**: [Developer name or AI assistant]

## Environment Information
**Operating System**: [e.g., Ubuntu 22.04, macOS 13.0, Windows 11]
**Python Version**: [e.g., 3.10.12]
**SARdine Version**: [e.g., v0.2.1]
**Dependencies**: [Relevant package versions]
**Hardware**: [CPU, Memory, GPU if relevant]

## Bug Description

### What Happened?
[Clear description of the observed behavior]

### What Should Have Happened?
[Description of expected behavior]

### Impact Assessment
- **Users Affected**: [Who is impacted]
- **Functionality Impact**: [What features are broken]
- **Data Impact**: [Is data integrity affected]
- **Performance Impact**: [Speed/memory degradation]

## Reproduction Steps
1. [First step to reproduce]
2. [Second step to reproduce]
3. [Continue with detailed steps]
4. [Final step that triggers the bug]

### Minimal Reproduction Case
```python
# Provide minimal code that reproduces the issue
# Include necessary imports and setup
import sardine

# Minimal example that triggers the bug
result = sardine.process_data(problematic_input)
```

### Test Data Requirements
- **Input Data**: [Description of data needed to reproduce]
- **Data Size**: [File sizes, array dimensions]
- **Data Location**: [Where to find test data]
- **Special Conditions**: [Specific data characteristics needed]

## Error Information

### Error Messages
```
[Copy and paste exact error messages here]
```

### Stack Trace
```
[Full stack trace if available]
```

### Log Output
```
[Relevant log entries leading up to the error]
```

## Investigation Notes

### Root Cause Analysis
[Preliminary analysis of what might be causing the issue]

### Code Location
- **File**: [Specific file where bug likely occurs]
- **Function/Method**: [Specific function]
- **Line Numbers**: [Approximate line numbers]

### Related Code Areas
- [Other files/functions that might be related]
- [Dependencies that could be involved]

## Diagnostic Information

### System State
- **Memory Usage**: [At time of error]
- **CPU Usage**: [At time of error]
- **Disk Space**: [Available space]
- **Network Status**: [If network-related]

### Configuration
```yaml
# Relevant configuration files or settings
# that might affect the bug
```

## Workaround
[If a temporary workaround exists, describe it here]

### Workaround Code
```python
# If applicable, provide code for temporary workaround
```

## Fix Strategy

### Proposed Solution
[Initial thoughts on how to fix the issue]

### Implementation Plan
1. [ ] Step 1: [Specific implementation step]
2. [ ] Step 2: [Another implementation step]
3. [ ] Step 3: [Testing and validation step]

### Alternative Approaches
- **Approach 1**: [Alternative solution]
  - **Pros**: [Benefits]
  - **Cons**: [Drawbacks]
- **Approach 2**: [Another alternative]
  - **Pros**: [Benefits]
  - **Cons**: [Drawbacks]

## Testing Strategy

### Unit Tests
- [ ] Test case 1: [Specific test for the fix]
- [ ] Test case 2: [Edge case test]
- [ ] Test case 3: [Regression test]

### Integration Tests
- [ ] Integration scenario 1: [Full pipeline test]
- [ ] Integration scenario 2: [Cross-module test]

### Regression Testing
- [ ] Verify existing functionality still works
- [ ] Test with various input data types
- [ ] Performance regression testing

## Validation Criteria
- [ ] Bug is completely resolved
- [ ] No new bugs introduced
- [ ] Performance is not degraded
- [ ] All existing tests pass
- [ ] New tests added for the fix
- [ ] Documentation updated if needed

## Related Issues
- **Similar Issues**: [Links to related bug reports]
- **Dependencies**: [Other issues that must be fixed first]
- **Blocks**: [Other issues that this bug is blocking]

## Communication Log

### Initial Report
- **Reporter**: [Who reported the bug]
- **Date**: [When it was reported]
- **Context**: [How it was discovered]

### Investigation Updates
- **Date**: [Update date]
  - **Finding**: [What was discovered]
  - **Next Steps**: [What to do next]

### Resolution
- **Fix Applied**: [Date and description of fix]
- **Verification**: [How the fix was verified]
- **Deployment**: [When and how fix was deployed]

## Lessons Learned
- **Prevention**: [How to prevent similar bugs]
- **Detection**: [How to catch similar issues earlier]
- **Process Improvement**: [Process changes needed]

---

**Reported by**: [Name]
**Date Reported**: [Report date]
**Last Updated**: [Last update date]
**Resolution Date**: [When fixed]
**Time to Resolution**: [Duration from report to fix]