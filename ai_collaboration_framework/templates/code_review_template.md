# Code Review Template

## Review Overview
**Pull Request/Commit**: [PR number or commit hash]
**Author**: [Code author name]
**Reviewer**: [Reviewer name - human or AI assistant]
**Review Date**: [Date of review]
**Review Type**: [Feature/Bug Fix/Refactor/Documentation/Other]

## Code Quality Assessment

### Functionality ✅❌
- [ ] **Correctness**: Code does what it's supposed to do
- [ ] **Completeness**: All requirements are implemented
- [ ] **Edge Cases**: Handles boundary conditions properly
- [ ] **Error Handling**: Appropriate exception handling
- [ ] **Input Validation**: Validates inputs correctly

**Score**: [1-5] | **Comments**: [Specific feedback]

### Code Structure ✅❌
- [ ] **Organization**: Logical file and function organization
- [ ] **Modularity**: Well-separated concerns
- [ ] **Reusability**: Code can be reused appropriately
- [ ] **Coupling**: Minimal unnecessary dependencies
- [ ] **Cohesion**: Related functionality grouped together

**Score**: [1-5] | **Comments**: [Specific feedback]

### Readability ✅❌
- [ ] **Naming**: Clear, descriptive variable/function names
- [ ] **Comments**: Adequate and helpful comments
- [ ] **Documentation**: Docstrings and inline docs
- [ ] **Formatting**: Consistent code formatting
- [ ] **Complexity**: Code complexity is manageable

**Score**: [1-5] | **Comments**: [Specific feedback]

### Performance ✅❌
- [ ] **Efficiency**: Algorithms are appropriately efficient
- [ ] **Memory Usage**: Reasonable memory consumption
- [ ] **Scalability**: Handles expected data volumes
- [ ] **Optimization**: No obvious performance bottlenecks
- [ ] **Resource Management**: Proper cleanup of resources

**Score**: [1-5] | **Comments**: [Specific feedback]

### Testing ✅❌
- [ ] **Unit Tests**: Adequate unit test coverage
- [ ] **Integration Tests**: Integration scenarios covered
- [ ] **Test Quality**: Tests are well-written and meaningful
- [ ] **Edge Cases**: Tests cover boundary conditions
- [ ] **Mock Usage**: Appropriate use of mocks/stubs

**Score**: [1-5] | **Comments**: [Specific feedback]

## Scientific/Technical Accuracy

### Algorithm Implementation ✅❌
- [ ] **Mathematical Correctness**: Formulas implemented correctly
- [ ] **Scientific Accuracy**: Follows established scientific principles
- [ ] **Literature Compliance**: Matches referenced papers/standards
- [ ] **Numerical Stability**: Handles floating-point operations safely
- [ ] **Domain Knowledge**: Demonstrates understanding of SAR processing

**Score**: [1-5] | **Comments**: [Specific feedback]

### Data Handling ✅❌
- [ ] **Data Integrity**: Preserves data throughout processing
- [ ] **Format Compliance**: Follows expected data formats
- [ ] **Metadata Handling**: Properly manages metadata
- [ ] **Coordinate Systems**: Correct spatial reference handling
- [ ] **Quality Flags**: Appropriate quality indicators

**Score**: [1-5] | **Comments**: [Specific feedback]

## Security and Safety

### Security Considerations ✅❌
- [ ] **Input Sanitization**: Protects against malicious input
- [ ] **File Operations**: Safe file handling practices
- [ ] **Dependencies**: Secure use of external libraries
- [ ] **Sensitive Data**: Proper handling of sensitive information
- [ ] **Access Control**: Appropriate permission handling

**Score**: [1-5] | **Comments**: [Specific feedback]

### Safety and Robustness ✅❌
- [ ] **Error Recovery**: Graceful handling of failures
- [ ] **Resource Limits**: Prevents resource exhaustion
- [ ] **Backwards Compatibility**: Maintains API compatibility
- [ ] **Deprecation**: Proper deprecation of old features
- [ ] **Logging**: Adequate logging for debugging

**Score**: [1-5] | **Comments**: [Specific feedback]

## Detailed Code Review

### Strengths
1. **[Strength 1]**: [Specific positive aspect]
2. **[Strength 2]**: [Another positive aspect]
3. **[Strength 3]**: [Additional strength]

### Areas for Improvement
1. **[Issue 1]**: [Specific issue description]
   - **Location**: [File:line or function name]
   - **Suggestion**: [How to improve]
   - **Priority**: [High/Medium/Low]

2. **[Issue 2]**: [Another issue description]
   - **Location**: [File:line or function name]
   - **Suggestion**: [How to improve]
   - **Priority**: [High/Medium/Low]

### Critical Issues (Must Fix)
- **Issue**: [Critical problem description]
  - **Impact**: [Why this is critical]
  - **Solution**: [Required fix]

### Suggestions (Nice to Have)
- **Suggestion**: [Improvement suggestion]
  - **Benefit**: [Why this would be helpful]
  - **Implementation**: [How to implement]

## Specific Line Comments

### File: [filename]
```python
# Line XX-YY: [Comment about specific code section]
def example_function():
    # Your specific feedback about this code
    pass
```

**Feedback**: [Detailed feedback about this section]

### File: [another_filename]
```python
# Line XX-YY: [Comment about another code section]
class ExampleClass:
    # Your specific feedback about this code
    pass
```

**Feedback**: [Detailed feedback about this section]

## Testing Validation

### Test Coverage Analysis
- **Current Coverage**: [Percentage]
- **Missing Coverage**: [Areas not covered]
- **Test Quality**: [Assessment of test effectiveness]

### Recommended Additional Tests
1. **Test Case**: [Description of needed test]
   - **Purpose**: [Why this test is needed]
   - **Implementation**: [How to implement]

2. **Test Case**: [Another needed test]
   - **Purpose**: [Why this test is needed]
   - **Implementation**: [How to implement]

## Documentation Review

### Code Documentation ✅❌
- [ ] **Docstrings**: Complete and accurate function/class docs
- [ ] **Type Hints**: Appropriate type annotations
- [ ] **Comments**: Helpful inline comments
- [ ] **README Updates**: Documentation reflects changes
- [ ] **API Documentation**: Public interfaces documented

### Documentation Quality
- **Clarity**: [Assessment of documentation clarity]
- **Completeness**: [Coverage of all functionality]
- **Examples**: [Quality and completeness of examples]

## Performance Analysis

### Performance Testing
- **Benchmarks Run**: [What performance tests were executed]
- **Results**: [Performance measurements]
- **Comparison**: [How it compares to previous versions]

### Performance Recommendations
- **Optimization 1**: [Specific performance improvement]
- **Optimization 2**: [Another performance suggestion]

## Integration Compatibility

### API Changes
- **Breaking Changes**: [Any breaking changes introduced]
- **New Features**: [New API features added]
- **Deprecations**: [Any features marked as deprecated]

### Dependency Impact
- **New Dependencies**: [Any new dependencies added]
- **Version Updates**: [Dependency version changes]
- **Compatibility**: [Impact on existing integrations]

## Overall Assessment

### Summary Score: [1-5]
**Overall Quality**: [Brief assessment]

### Recommendation
- [ ] **Approve**: Ready to merge as-is
- [ ] **Approve with Minor Changes**: Can merge after small fixes
- [ ] **Request Changes**: Needs significant improvements
- [ ] **Reject**: Major issues require complete rework

### Next Steps
1. **Immediate Actions**: [What needs to be done right away]
2. **Follow-up Tasks**: [Future improvements to consider]
3. **Monitoring**: [What to watch after deployment]

## Collaboration Notes

### Learning Points
- **For Author**: [Key learning points for code author]
- **For Project**: [Insights for future development]
- **Best Practices**: [Patterns to adopt going forward]

### Knowledge Transfer
- **Technical Insights**: [Important technical knowledge shared]
- **Domain Knowledge**: [Scientific/domain insights]
- **Process Improvements**: [Process lessons learned]

---

**Review Completed by**: [Reviewer name]
**Review Duration**: [Time spent on review]
**Follow-up Required**: [Yes/No - describe if yes]
**Approval Status**: [Approved/Changes Requested/Rejected]