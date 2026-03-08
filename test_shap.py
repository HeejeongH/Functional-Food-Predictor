import requests
import json

print("="*60)
print("🧪 Testing SHAP Analysis with Improved Service")
print("="*60)

url = "http://localhost:3000/api/shap/analyze"
payload = {
    "model_id": "FTO_MLModelType.XGBOOST_20260306_053810",
    "feature_type": "fingerprint",
    "top_n": 10
}

print(f"\n📤 Sending POST request to {url}")
print(f"   Payload: {json.dumps(payload, indent=2)}")

try:
    response = requests.post(url, json=payload, timeout=600)  # 10분 타임아웃
    
    print(f"\n📥 Response Status: {response.status_code}")
    
    if response.status_code == 200:
        result = response.json()
        print(f"\n✅ SHAP Analysis Successful!")
        print(f"\n🏆 Top 10 Features:")
        for i, feat in enumerate(result['top_features'][:10], 1):
            print(f"   {i:2d}. {feat['feature']}: {feat['importance']:.6f}")
        
        print(f"\n📊 Summary:")
        summary = result['shap_values_summary']
        print(f"   Mean |SHAP|: {summary['mean_abs_shap']:.6f}")
        print(f"   Max |SHAP|:  {summary['max_abs_shap']:.6f}")
        print(f"   Samples:     {summary['samples_analyzed']}")
        print(f"\n🖼️  Plot: {result['plot_path']}")
    else:
        print(f"\n❌ Error {response.status_code}:")
        print(f"   {response.text}")

except requests.exceptions.Timeout:
    print("\n❌ Request timed out after 10 minutes")
except Exception as e:
    print(f"\n❌ Error: {str(e)}")

print(f"\n{'='*60}\n")
