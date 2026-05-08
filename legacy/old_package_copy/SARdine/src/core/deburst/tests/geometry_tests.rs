use crate::core::deburst::geometry::{
    build_line_timing_with_offset, eval_dc_fm, extract_dc_estimates_from_annotation,
    extract_fm_estimates_from_annotation, DcEstimate, RangePolynomial,
};

#[test]
fn line_timing_with_offset_is_centered() {
    let lines = 10;
    let az_interval = 0.001;
    let time_offset = 5.0;

    let timings = build_line_timing_with_offset(lines, az_interval, time_offset);
    assert_eq!(timings.len(), lines);

    let center = (lines as f64 - 1.0) * 0.5;
    let expected_first = (0.0 - center) * az_interval + time_offset;
    assert!((timings[0].t_az - expected_first).abs() < 1e-9);
    let center_idx = lines / 2;
    assert!((timings[center_idx].t_az - time_offset).abs() < az_interval);
}

#[test]
fn dc_and_fm_evaluation_matches_coeffs() {
    let dc_poly = vec![100.0, 0.0];
    let fm_poly = vec![2000.0, 0.0];
    let (dc, fm) = eval_dc_fm(0.001, &dc_poly, &fm_poly, None, None);
    assert!((dc - 100.0).abs() < 1e-6);
    assert!((fm - 2000.0).abs() < 1e-6);
}

#[test]
fn extract_dc_estimates_with_t0() {
    let xml = r#"
    <dopplerCentroid>
        <dcEstimateList count=\"2\">
            <dcEstimate>
                <t0>0.5</t0>
                <slantRangeTime>5.5e-3</slantRangeTime>
                <dataDcPolynomial>1.0 2.0 3.0</dataDcPolynomial>
            </dcEstimate>
            <dcEstimate>
                <t0>-1.25</t0>
                <slantRangeTime>5.9e-3</slantRangeTime>
                <dataDcPolynomial>4 5 6 7</dataDcPolynomial>
            </dcEstimate>
        </dcEstimateList>
    </dopplerCentroid>
    "#;

    let estimates = extract_dc_estimates_from_annotation(xml).expect("should parse dc estimates");

    assert_eq!(estimates.len(), 2);
    assert_eq!(estimates[0].coeffs, vec![1.0, 2.0, 3.0]);
    assert_eq!(estimates[0].t0, Some(0.5));
    assert_eq!(estimates[0].slant_range_time, Some(5.5e-3));
    assert_eq!(estimates[1].coeffs, vec![4.0, 5.0, 6.0, 7.0]);
    assert_eq!(estimates[1].t0, Some(-1.25));
    assert_eq!(estimates[1].slant_range_time, Some(5.9e-3));
}

#[test]
fn extract_fm_estimates_with_t0() {
    let xml = r#"
    <fmRateList count=\"1\">
        <fmRate>
            <t0>2.5</t0>
            <slantRangeTime>5.7e-3</slantRangeTime>
            <dataFmratePolynomial>10 20 30</dataFmratePolynomial>
        </fmRate>
    </fmRateList>
    "#;

    let estimates = extract_fm_estimates_from_annotation(xml).expect("should parse fm estimates");

    assert_eq!(estimates.len(), 1);
    assert_eq!(estimates[0].coeffs, vec![10.0, 20.0, 30.0]);
    assert_eq!(estimates[0].t0, Some(2.5));
    assert_eq!(estimates[0].slant_range_time, Some(5.7e-3));
}

#[test]
fn range_polynomial_interpolation() {
    let estimates = vec![
        DcEstimate {
            coeffs: vec![1.0, 0.0],
            t0: None,
            azimuth_time: None,
            slant_range_time: Some(5.0e-3),
        },
        DcEstimate {
            coeffs: vec![3.0, 0.0],
            t0: None,
            azimuth_time: None,
            slant_range_time: Some(6.0e-3),
        },
    ];

    let model = RangePolynomial::from_dc_estimates(&estimates).expect("should build model");
    let v_mid = model.evaluate(5.5e-3, None);
    assert!((v_mid - 2.0).abs() < 1e-9);
    let v_low = model.evaluate(4.0e-3, None);
    assert!((v_low - 1.0).abs() < 1e-9);
    let v_high = model.evaluate(6.5e-3, None);
    assert!((v_high - 3.0).abs() < 1e-9);
}
