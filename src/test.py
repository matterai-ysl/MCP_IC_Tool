conv_analysis ={'force_convergence': {'converged': True, 'threshold': 0.01, 'final_max_force': 0.0, 'force_history': [0.0, 0.0, 0.0]}, 'energy_convergence': {'converged': False, 'final_energy': -17.5793482, 'energy_history': [25.60111643, -17.85724782, -17.90338016, -17.90340472, -17.90340472, -17.46948842, -17.54917201, -17.56080674, -17.57654541, -17.57812065, -17.57920307, -17.57929796, -17.57930134, -17.57930134, -17.57960382, -17.5793681, -17.57934883, -17.57934739, -17.57934739, -17.58025384, -17.57942086, -17.57935399, -17.5793482, -17.5793482], 'energy_changes': [43.45836425, 0.04613233999999977, 2.4559999999951287e-05, 0.0, 0.4339162999999999, 0.07968358999999836, 0.011634730000000815, 0.015738670000001065, 0.0015752399999975353, 0.0010824199999994732, 9.489000000328929e-05, 3.3799999989980734e-06, 0.0, 0.00030247999999843955, 0.00023571999999916216, 1.9269999999238507e-05, 1.4400000019065828e-06, 0.0, 0.0009064500000022235, 0.0008329800000019816, 6.68699999977207e-05, 5.790000003003115e-06, 0.0], 'avg_recent_change': 0.0003624180000009858}, 'electronic_convergence': {}, 'tail_check': {'matched': True, 'keywords': ['reached required accuracy', 'Voluntary'], 'tail_bytes': 1024, 'exception': None}, 'overall_convergence': True}



print(conv_analysis['force_convergence']['converged'])
print(conv_analysis['energy_convergence']['converged'])
print(conv_analysis['electronic_convergence'])
print(conv_analysis['tail_check'])

result = {
    'force_convergence': conv_analysis.get('force_convergence', {}).get("converged", False),
    'final_max_force': conv_analysis.get('force_convergence', {}).get("final_max_force", None),
    'energy_convergence': conv_analysis.get('energy_convergence', {}).get("converged", False),
    'final_energy': conv_analysis.get('energy_convergence', {}).get("final_energy", None)
}
print("--------------------------------")
print(result)
print("--------------------------------")