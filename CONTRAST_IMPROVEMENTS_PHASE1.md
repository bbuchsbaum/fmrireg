# Phase 1 Improvements to contrast.R

## Summary of Changes

### 1. Added Internal Constants
- Created `.CONTRAST_TOLERANCE <- 1e-8` to replace magic numbers
- Updated all hardcoded tolerance values to use the constant

### 2. Consolidated Duplicate Functions
- Merged `.mask_to_weights()` and `.make_weights()` into single `.calculate_mask_weights()`
- Added backward compatibility aliases to prevent breaking existing code
- Improved error messages to use consistent function name

### 3. Removed Dead Code
- Removed unused `.col_index()` function
- Removed commented-out code on line 18

### 4. Enhanced Input Validation
- Added comprehensive validation to all constructor functions:
  - `contrast()`: validates form, name, and where parameters
  - `unit_contrast()`: validates A, name, and where parameters  
  - `pair_contrast()`: validates A, B, name, and where parameters
  - `oneway_contrast()`: validates A, name, and where parameters
  - `interaction_contrast()`: validates A, name, and where parameters
  - `poly_contrast()`: validates A, name, where, degree, and value_map parameters
  - `column_contrast()`: validates pattern_A, pattern_B, name parameters
- All validation includes helpful error messages

### 5. Improved Documentation
- Fixed `@export` directive on `translate_legacy_pattern()` (removed since it's internal)
- Added examples to `translate_legacy_pattern()`
- Enhanced examples for `contrast()`, `unit_contrast()`, and `pair_contrast()`

### 6. Code Quality Improvements
- Consistent error message formatting
- Better parameter validation with descriptive messages
- Improved code organization

## Benefits
- **No Breaking Changes**: All improvements maintain backward compatibility
- **Better Error Messages**: Users get clearer feedback when providing invalid inputs
- **Reduced Duplication**: Easier to maintain with consolidated weight calculation
- **Improved Testability**: Constants and validation make testing more reliable
- **Enhanced Usability**: Better documentation and examples help users

## Next Steps (Phase 2-5)
1. **Phase 2**: Reorganize code - move internal functions to separate file
2. **Phase 3**: Performance improvements - cache regex patterns, vectorize operations
3. **Phase 4**: Add validation methods and utility functions
4. **Phase 5**: Comprehensive documentation with more examples

## Testing
All existing tests continue to pass, confirming backward compatibility is maintained.